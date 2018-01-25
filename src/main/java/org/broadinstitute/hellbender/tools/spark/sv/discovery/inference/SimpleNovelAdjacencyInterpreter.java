package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple3;

import java.util.*;
import java.util.stream.Collectors;

/**
 * This deals with the special case where a contig has exactly two alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromMultipleAlignments()}.
 *
 * TODO: 1/19/18 see ticket 4189
 *      Exactly how the returned type in {@link SimpleNovelAdjacencyAndChimericAlignmentEvidence} is treated (trusted, or updated, or re-interpreted),
 *      is to be developed.
 */
public final class SimpleNovelAdjacencyInterpreter {

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public JavaPairRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>> inferTypeFromSingleContigSimpleChimera(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                                                                                              final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;

        final JavaRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence> simpleNovelAdjacencies =
                getSimpleNovelAdjacencyJavaRDD(assemblyContigs, svDiscoveryInputData);

        return simpleNovelAdjacencies
                        .mapToPair(simpleNovelAdjacencyAndChimericAlignmentEvidence ->
                                new Tuple2<>(simpleNovelAdjacencyAndChimericAlignmentEvidence,
                                        inferSimpleOrBNDTypesFromNovelAdjacency(simpleNovelAdjacencyAndChimericAlignmentEvidence,
                                                referenceBroadcast.getValue(), referenceSequenceDictionaryBroadcast.getValue())));
    }

    private JavaRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence> getSimpleNovelAdjacencyJavaRDD(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                                                                     final SvDiscoveryInputData svDiscoveryInputData) {
        final Logger toolLogger = svDiscoveryInputData.toolLogger;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputData.discoverStageArgs;
        final List<SVInterval> assembledIntervals = svDiscoveryInputData.assembledIntervals;

        final JavaRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence> simpleNovelAdjacencies =
                assemblyContigs
                        .filter(tig ->
                                ChimericAlignment.splitPairStrongEnoughEvidenceForCA(tig.getSourceContig().alignmentIntervals.get(0),
                                        tig.getSourceContig().alignmentIntervals.get(1),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ, MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(tig -> {
                            final SAMSequenceDictionary refSeqDict = referenceSequenceDictionaryBroadcast.getValue();
                            final ChimericAlignment simpleChimera = ChimericAlignment.extractSimpleChimera(tig, refSeqDict);
                            final byte[] contigSequence = tig.getSourceContig().contigSequence;

                            final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype =
                                    new NovelAdjacencyAndInferredAltHaptype(simpleChimera, contigSequence, refSeqDict);
                            return new Tuple2<>(novelAdjacencyAndInferredAltHaptype, simpleChimera);
                        })
                        .groupByKey()       // group the same novel adjacency produced by different contigs together
                        .map(noveltyAndEvidence ->
                                new SimpleNovelAdjacencyAndChimericAlignmentEvidence(noveltyAndEvidence._1, Lists.newArrayList(noveltyAndEvidence._2)));

        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals,
                simpleNovelAdjacencies.map(SimpleNovelAdjacencyAndChimericAlignmentEvidence::getNovelAdjacencyReferenceLocations).collect(),
                referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);

        return simpleNovelAdjacencies;
    }

    /**
     * Main method: given simple novel adjacency
     *  (that is, affected reference locations, alt haplotype sequence, and chimeric alignment evidence),
     *  infer type.
     *
     * @return the inferred type could be a single entry for simple variants, or a list of two entries with BND mates.
     */
    static List<SvType> inferSimpleOrBNDTypesFromNovelAdjacency(final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence,
                                                                final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary) {

        // based on characteristic of simple chimera, infer type
        final List<ChimericAlignment> alignmentEvidence = simpleNovelAdjacencyAndChimericAlignmentEvidence.getAlignmentEvidence();

        final List<SvType> inferredType;
        final NovelAdjacencyAndInferredAltHaptype novelAdjacency = simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations();
        if ( alignmentEvidence.stream().allMatch(ChimericAlignment::isLikelySimpleTranslocation) ) { // all indicate simple translocation
            final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForTranslocSuspect =
                    BreakEndVariantType.TransLocBND.getOrderedMates(novelAdjacency,
                            reference, referenceDictionary);
            inferredType = Arrays.asList(orderedMatesForTranslocSuspect._1, orderedMatesForTranslocSuspect._2);
        } else if ( alignmentEvidence.stream().allMatch(ChimericAlignment::isLikelyInvertedDuplication) ) { // all indicate inverted duplication
            inferredType = Collections.singletonList( new SimpleSVType.DuplicationInverted(novelAdjacency) );
        } else if ( alignmentEvidence.stream().map(ca -> ca.strandSwitch).noneMatch(ss -> ss.equals(StrandSwitch.NO_SWITCH)) ) { // all indicate simple (i.e. no duplicate) strand-switch novel adjacency
            final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForInversionSuspect =
                    BreakEndVariantType.InvSuspectBND.getOrderedMates(novelAdjacency, reference);
            inferredType = Arrays.asList(orderedMatesForInversionSuspect._1, orderedMatesForInversionSuspect._2);
        } else if ( alignmentEvidence.stream().allMatch(ChimericAlignment::isNeitherSimpleTranslocationNorIncompletePicture) &&
                alignmentEvidence.stream().map(ca -> ca.strandSwitch).allMatch(ss -> ss.equals(StrandSwitch.NO_SWITCH)) ){ // all point to simple insertion/deletion/small duplication
            inferredType = Collections.singletonList( inferSimpleTypeFromNovelAdjacency(novelAdjacency) );
        } else {
            throw new GATKException
                    .ShouldNeverReachHereException("novel adjacency has its supporting chimeric alignments showing inconsistent behavior\n" +
                    simpleNovelAdjacencyAndChimericAlignmentEvidence.toString());
        }

        return inferredType;
    }

    @VisibleForTesting
    public static SimpleSVType inferSimpleTypeFromNovelAdjacency(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {

        final int start = novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd();
        final int end = novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart();
        final StrandSwitch strandSwitch = novelAdjacencyAndInferredAltHaptype.strandSwitch;

        final SimpleSVType type;
        if (strandSwitch == StrandSwitch.NO_SWITCH) { // no strand switch happening, so no inversion
            if (start==end) { // something is inserted
                final boolean hasNoDupSeq = !novelAdjacencyAndInferredAltHaptype.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyAndInferredAltHaptype.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        throw new GATKException("Something went wrong in type inference, there's suspected insertion happening but no inserted sequence could be inferred "
                                + novelAdjacencyAndInferredAltHaptype.toString());
                    } else {
                        type = new SimpleSVType.Insertion(novelAdjacencyAndInferredAltHaptype); // simple insertion (no duplication)
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.DuplicationTandem(novelAdjacencyAndInferredAltHaptype); // clean expansion of repeat 1 -> 2, or complex expansion
                    } else {
                        type = new SimpleSVType.DuplicationTandem(novelAdjacencyAndInferredAltHaptype); // expansion of 1 repeat on ref to 2 repeats on alt with inserted sequence in between the 2 repeats
                    }
                }
            } else {
                final boolean hasNoDupSeq = !novelAdjacencyAndInferredAltHaptype.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyAndInferredAltHaptype.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyAndInferredAltHaptype); // clean deletion
                    } else {
                        type = new SimpleSVType.Deletion(novelAdjacencyAndInferredAltHaptype); // scarred deletion
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyAndInferredAltHaptype); // clean contraction of repeat 2 -> 1, or complex contraction
                    } else {
                        throw new GATKException("Something went wrong in novel adjacency interpretation: " +
                                " inferring simple SV type from a novel adjacency between two different reference locations, but annotated with both inserted sequence and duplication, which is NOT simple.\n"
                                + novelAdjacencyAndInferredAltHaptype.toString());
                    }
                }
            }
        } else {
            type = new SimpleSVType.Inversion(novelAdjacencyAndInferredAltHaptype);
        }

        return type;
    }

    /**
     * Infers exact position of breakpoints based on characteristics of input {@link ChimericAlignment},
     * and alt haplotype sequence based on given contig sequence.
     *
     * <p>
     *     If there is homologous sequence represented in the {@link ChimericAlignment}, it will be assigned to
     *     the side of the breakpoint with higher reference coordinates (as judged by {@link SAMSequenceDictionary}),
     *     i.e. we follow left alignment convention.
     * </p>
     */
    abstract static class BreakpointsInference {

        protected String upstreamBreakpointRefContig;
        protected String downstreamBreakpointRefContig;
        protected int upstreamBreakpointRefPos;
        protected int downstreamBreakpointRefPos;

        final Tuple2<SimpleInterval, SimpleInterval> getLeftJustifiedBreakpoints() {
            return new Tuple2<>(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                                new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos));
        }

        // TODO: 1/24/18 this is the task of next commit
        final byte[] getInferredAltHaplotypeSequence() {
            return null;
        }

        static final BreakpointsInference getInferenceClass(final ChimericAlignment chimericAlignment,
                                                            final BreakpointComplications complication,
                                                            final SAMSequenceDictionary referenceDictionary) {
            final boolean sameChromosome =
                    chimericAlignment.regionWithLowerCoordOnContig.referenceSpan.getContig()
                            .equals(chimericAlignment.regionWithHigherCoordOnContig.referenceSpan.getContig());
            if (sameChromosome) {
                if (complication.hasDuplicationAnnotation()) {
                    if (complication.hasDupSeqButNoStrandSwitch()) { // tandem duplication
                        return new TanDupBreakpointsInference(chimericAlignment, complication);
                    } else { // inverted duplication
                        return new InvDupBreakpointsInference(chimericAlignment, complication);
                    }
                } else { // no explicit dup annotation (ins/del, dispersed dup, inv)
                    return new NoExplicitDupAnnotationBreakpointsInference(chimericAlignment, complication);
                }
            } else {
                return new TransLocBreakpointsAligner(chimericAlignment, complication, referenceDictionary);
            }
        }

        ///////////////
        abstract static class SameChrEventsBreakpointsInference extends BreakpointsInference {
            protected SameChrEventsBreakpointsInference(final ChimericAlignment ca) {
                upstreamBreakpointRefContig = downstreamBreakpointRefContig
                        = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
            }
        }

        /**
         * For simple insertion/deletion, inversion, and (suspected) dispersed duplications without explicit duplication annotation
         */
        final static class NoExplicitDupAnnotationBreakpointsInference extends SameChrEventsBreakpointsInference {
            NoExplicitDupAnnotationBreakpointsInference(final ChimericAlignment ca, final BreakpointComplications complication) {

                super(ca);

                final int homologyLen = complication.getHomologyForwardStrandRep().length();
                final AlignmentInterval one = ca.regionWithLowerCoordOnContig,
                                        two = ca.regionWithHigherCoordOnContig;
                final SimpleInterval leftReferenceInterval, rightReferenceInterval;
                if (ca.isForwardStrandRepresentation) {
                    leftReferenceInterval  = one.referenceSpan;
                    rightReferenceInterval = two.referenceSpan;
                } else {
                    leftReferenceInterval  = two.referenceSpan;
                    rightReferenceInterval = one.referenceSpan;
                }
                if (ca.strandSwitch == StrandSwitch.NO_SWITCH) {

                    final List<AlignmentInterval> deOverlappedTempAlignments = ContigAlignmentsModifier.removeOverlap(one, two, null);
                    final AlignmentInterval deOverlappedOne = deOverlappedTempAlignments.get(0),
                                            deOverlappedTwo = deOverlappedTempAlignments.get(1);
                    final boolean isLikelyDispersedDuplication =
                            deOverlappedOne.referenceSpan.getStart() > deOverlappedTwo.referenceSpan.getStart() == deOverlappedOne.forwardStrand;
                    if (isLikelyDispersedDuplication) { // left and right seem to be flipped but that's the feature of dispersed duplication
                        upstreamBreakpointRefPos = rightReferenceInterval.getStart();
                        downstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                    } else {
                        upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                        downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;
                    }

                } else if (ca.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE){
                    upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                    downstreamBreakpointRefPos = rightReferenceInterval.getEnd();
                } else {
                    upstreamBreakpointRefPos = leftReferenceInterval.getStart() - 1;
                    downstreamBreakpointRefPos = rightReferenceInterval.getStart() + homologyLen - 1;
                }

            }
        }

        final static class TanDupBreakpointsInference extends SameChrEventsBreakpointsInference {

            TanDupBreakpointsInference(final ChimericAlignment ca,
                                       final BreakpointComplications complication) {
                super(ca);

                final int homologyLen = complication.getHomologyForwardStrandRep().length();

                final SimpleInterval leftReferenceInterval, rightReferenceInterval;
                if (ca.isForwardStrandRepresentation) {
                    leftReferenceInterval  = ca.regionWithLowerCoordOnContig.referenceSpan;
                    rightReferenceInterval = ca.regionWithHigherCoordOnContig.referenceSpan;
                } else {
                    leftReferenceInterval  = ca.regionWithHigherCoordOnContig.referenceSpan;
                    rightReferenceInterval = ca.regionWithLowerCoordOnContig.referenceSpan;
                }
                if (complication.getDupSeqRepeatNumOnCtg() > complication.getDupSeqRepeatNumOnRef()) {
                    upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen
                            - (complication.getDupSeqRepeatNumOnCtg() - complication.getDupSeqRepeatNumOnRef()) * complication.getDupSeqRepeatUnitRefSpan().size();
                } else {
                    upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                }
                downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;;
            }
        }

        final static class InvDupBreakpointsInference extends SameChrEventsBreakpointsInference {
            InvDupBreakpointsInference(final ChimericAlignment ca,
                                       final BreakpointComplications complication) {
                super(ca);
                upstreamBreakpointRefPos = complication.getDupSeqRepeatUnitRefSpan().getStart() - 1;
                downstreamBreakpointRefPos = complication.getDupSeqRepeatUnitRefSpan().getEnd();
            }

            // TODO: 1/21/18 hookup at the right place (right now no variants are using this any way because inverted duplication contigs are filtered out)
            static Iterator<Tuple2<Tuple3<NovelAdjacencyAndInferredAltHaptype, SimpleSVType.DuplicationInverted, byte[]>, List<ChimericAlignment>>>
            inferInvDupRange(final Tuple2<NovelAdjacencyAndInferredAltHaptype, Iterable<Tuple2<ChimericAlignment, byte[]>>> noveltyAndEvidence) {

                final NovelAdjacencyAndInferredAltHaptype novelAdjacency = noveltyAndEvidence._1;
                final SimpleSVType.DuplicationInverted duplicationInverted = new SimpleSVType.DuplicationInverted(novelAdjacency);

                // doing this because the same novel adjacency reference locations might be induced by different (probably only slightly) alt haplotypes, so a single group by NARL is not enough
                final Iterable<Tuple2<ChimericAlignment, byte[]>> chimeraAndContigSeq = noveltyAndEvidence._2;
                final Set<Map.Entry<byte[], List<ChimericAlignment>>> alignmentEvidenceGroupedByAltHaplotypeSequence =
                        Utils.stream(chimeraAndContigSeq)
                                .collect(
                                        Collectors.groupingBy(caAndSeq ->
                                                        extractAltHaplotypeForInvDup(caAndSeq._1, caAndSeq._2),
                                                Collectors.mapping(caAndSeq -> caAndSeq._1, Collectors.toList())
                                        )
                                )
                                .entrySet();

                return alignmentEvidenceGroupedByAltHaplotypeSequence.stream()
                        .map(entry -> new Tuple2<>(new Tuple3<>(novelAdjacency, duplicationInverted, entry.getKey()),
                                entry.getValue()))
                        .iterator();
            }

            /**
             * Extract alt haplotype sequence, based on the input {@code chimericAlignment} and {@code contigSeq}.
             */
            private static byte[] extractAltHaplotypeForInvDup(final ChimericAlignment chimericAlignment, final byte[] contigSeq) {

                final AlignmentInterval firstAlignmentInterval  = chimericAlignment.regionWithLowerCoordOnContig;
                final AlignmentInterval secondAlignmentInterval = chimericAlignment.regionWithHigherCoordOnContig;

                final int start, end; // intended to be 0-based, semi-open [start, end)
                final boolean needRC;
                // below we need to use cigars of the provided alignments to compute how long we need to walk on the read
                // so that we can "start" to or "end" to collect bases for alternative haplotype sequence,
                // because one could imagine either alignment has long flanking region that is far from the affected reference region.
                if (firstAlignmentInterval.forwardStrand) {
                    final int alpha = firstAlignmentInterval.referenceSpan.getStart(),
                            omega = secondAlignmentInterval.referenceSpan.getStart();
                    if (alpha <= omega) {
                        final int walkOnReadUntilDuplicatedSequence ;
                        if (alpha == omega) {
                            walkOnReadUntilDuplicatedSequence = 0;
                        } else {
                            walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(firstAlignmentInterval.cigarAlong5to3DirectionOfContig,
                                    firstAlignmentInterval.startInAssembledContig, omega - alpha, false);
                        }
                        start = firstAlignmentInterval.startInAssembledContig + walkOnReadUntilDuplicatedSequence - 1;
                        end = secondAlignmentInterval.endInAssembledContig;
                        needRC = false;
                    } else {
                        final int walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(secondAlignmentInterval.cigarAlong5to3DirectionOfContig,
                                secondAlignmentInterval.endInAssembledContig, alpha - omega, true);
                        start = firstAlignmentInterval.startInAssembledContig - 1;
                        end = secondAlignmentInterval.endInAssembledContig - walkOnReadUntilDuplicatedSequence;
                        needRC = true;
                    }
                } else {
                    final int alpha = firstAlignmentInterval.referenceSpan.getEnd(),
                            omega = secondAlignmentInterval.referenceSpan.getEnd();
                    if (alpha >= omega) {
                        final int walkOnReadUntilDuplicatedSequence ;
                        if (alpha == omega) {
                            walkOnReadUntilDuplicatedSequence = 0;
                        } else {
                            walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(firstAlignmentInterval.cigarAlong5to3DirectionOfContig,
                                    firstAlignmentInterval.startInAssembledContig, alpha - omega, false);
                        }
                        start = firstAlignmentInterval.startInAssembledContig + walkOnReadUntilDuplicatedSequence - 1;
                        end = secondAlignmentInterval.endInAssembledContig;
                        needRC = true;
                    } else {
                        final int walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(secondAlignmentInterval.cigarAlong5to3DirectionOfContig,
                                secondAlignmentInterval.endInAssembledContig, omega - alpha, true);
                        start = firstAlignmentInterval.startInAssembledContig - 1;
                        end = secondAlignmentInterval.endInAssembledContig - walkOnReadUntilDuplicatedSequence;
                        needRC = false;
                    }
                }

                final byte[] seq = Arrays.copyOfRange(contigSeq, start, end);
                if (needRC) SequenceUtil.reverseComplement(seq, 0, seq.length);
                return seq;
            }
        }

        ///////////////

        /**
         * For computing exact and left-adjusted breakpoint locations of inter-chromosome novel adjacency,
         * with or without strand switch.
         */
        final static class TransLocBreakpointsAligner extends BreakpointsInference {

            TransLocBreakpointsAligner(final ChimericAlignment ca,
                                       final BreakpointComplications complication,
                                       final SAMSequenceDictionary referenceDictionary) {


                determineRefContigs(ca, referenceDictionary);

                extractRefPositions(ca, complication, referenceDictionary);
            }

            private void extractRefPositions(final ChimericAlignment ca, final BreakpointComplications complication,
                                             final SAMSequenceDictionary referenceDictionary) {
                final int homologyLen = complication.getHomologyForwardStrandRep().length();
                final boolean firstInPartner = isFirstInPartner(ca, referenceDictionary);
                if (firstInPartner) {
                    switch (ca.strandSwitch) {
                        case NO_SWITCH:
                            if (ca.isForwardStrandRepresentation) {
                                upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                                downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                            } else {
                                upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                                downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            }
                            break;
                        case FORWARD_TO_REVERSE:
                            upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd();
                            break;
                        case REVERSE_TO_FORWARD:
                            upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                            downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart() + homologyLen;
                            break;
                        default: throw new GATKException("Unseen strand switch case for: " + ca.toString());
                    }
                } else {
                    switch (ca.strandSwitch) {
                        case NO_SWITCH:
                            if (ca.isForwardStrandRepresentation) {
                                upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                                downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            } else {
                                upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                                downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                            }
                            break;
                        case FORWARD_TO_REVERSE:
                            upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd();
                            break;
                        case REVERSE_TO_FORWARD:
                            upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                            downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart() + homologyLen;
                            break;
                        default: throw new GATKException("Unseen strand switch case for: " + ca.toString());
                    }
                }
            }

            private void determineRefContigs(ChimericAlignment ca, SAMSequenceDictionary referenceDictionary) {
                final boolean firstInPartner = isFirstInPartner(ca, referenceDictionary);
                if (firstInPartner) {
                    upstreamBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
                    downstreamBreakpointRefContig = ca.regionWithHigherCoordOnContig.referenceSpan.getContig();
                } else {
                    upstreamBreakpointRefContig = ca.regionWithHigherCoordOnContig.referenceSpan.getContig();
                    downstreamBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
                }
            }

            private static boolean isFirstInPartner(final ChimericAlignment ca, final SAMSequenceDictionary referenceDictionary) {
                switch (ca.strandSwitch) {
                    case NO_SWITCH: return 0 > IntervalUtils.compareContigs(ca.regionWithLowerCoordOnContig.referenceSpan,
                            ca.regionWithHigherCoordOnContig.referenceSpan, referenceDictionary);
                    case FORWARD_TO_REVERSE: case REVERSE_TO_FORWARD:
                        return ca.isForwardStrandRepresentation;
                    default:
                        throw new GATKException("Unseen strand switch case for: " + ca.toString());
                }
            }
        }
    }
}
