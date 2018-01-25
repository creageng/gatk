package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

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
                                new SimpleNovelAdjacencyAndChimericAlignmentEvidence(noveltyAndEvidence._1, noveltyAndEvidence._2));

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

}
