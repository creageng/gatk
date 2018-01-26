package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.Arrays;

/**
 * Based on alignment signature of the input simple chimera, infers
 * <ul>
 *     <li>
 *         exact position of breakpoints following the left-aligning convention,
 *     </li>
 *     <li>
 *         alt haplotype sequence based on given contig sequence
 *     </li>
 *     <li>
 *         complications such as homology, inserted sequence and duplicated ref region, if any.
 *         If there is homologous sequence represented in the {@link ChimericAlignment},
 *         it will be assigned to the side of the breakpoint with higher reference coordinates
 *         (as judged by {@link SAMSequenceDictionary}), i.e. we follow left alignment convention.
 *     </li>
 * </ul>
 */
abstract class BreakpointsInference {

    protected String upstreamBreakpointRefContig;
    protected String downstreamBreakpointRefContig;
    protected int upstreamBreakpointRefPos;
    protected int downstreamBreakpointRefPos;

    // TODO: 1/26/18 for alt haplotype sequence, we need a single policy to decide if the POS base should be included or not
    protected byte[] altHaplotypeSequence;
    protected BreakpointComplications complications;

    final Tuple2<SimpleInterval, SimpleInterval> getLeftJustifiedBreakpoints() {

        return
                new Tuple2<>(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                             new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos));
    }
    final byte[] getInferredAltHaplotypeSequence() {
        return altHaplotypeSequence;
    }
    final BreakpointComplications getComplications() {
        return complications;
    }

    ///////////////
    static final BreakpointsInference getInferenceClass(final ChimericAlignment simpleChimera,
                                                        final byte[] contigSequence,
                                                        final SAMSequenceDictionary referenceDictionary) {

        if ( simpleChimera.isLikelySimpleTranslocation() ) { // see {@link ChimericAlignment.isLikelySimpleTranslocation()} for definition
            final boolean sameChromosomeEvent =
                    simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig()
                            .equals(simpleChimera.regionWithHigherCoordOnContig.referenceSpan.getContig());
            if ( sameChromosomeEvent ) {
                return new IntraChrRefOrderSwapBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
            } else {
                return new InterChromosomeBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
            }
        } else {
            if (simpleChimera.strandSwitch != StrandSwitch.NO_SWITCH) {
                // TODO: 9/9/17 the case involves an inversion, could be retired once same chr strand-switch BND calls are evaluated.
                return new IntraChrStrandSwitchBreakpointInference(simpleChimera, contigSequence, referenceDictionary);
            } else {
                // note: here we are same chromosome, no ref order swap, no strand switch;
                // so for a simple chimera we are down to: ins, del, replacement, or small dup
                switch (inferFromSimpleChimera(simpleChimera)) {
                    case SIMPLE_DEL:
                        return new SimpleDeletionBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
                    case RPL:
                        return new ReplacementBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
                    case SIMPLE_INS:
                        return new InsertionBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
                    case SMALL_DUP_EXPANSION: case DEL_DUP_CONTRACTION: case SMALL_DUP_CPX:
                        return new SmallDuplicationBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
                    default:
                        throw new GATKException.ShouldNeverReachHereException("Inferred type not recognized for simple chimera:\t" + simpleChimera.toString());
                }
            }
        }
    }

    private enum TypeInferredFromSimpleChimera {
        SIMPLE_DEL, DEL_DUP_CONTRACTION, SIMPLE_INS, RPL, SMALL_DUP_EXPANSION, SMALL_DUP_CPX;
    }

    /**
     * Implementing a logic, where based on how the {@code simpleChimera} of a read/assembly contig overlap or
     * distant from each other, we infer possible simple type.
     * @throws IllegalArgumentException when the simple chimera indicates strand switch or simple translocation or incomplete picture.
     */
    static TypeInferredFromSimpleChimera inferFromSimpleChimera (final ChimericAlignment simpleChimera) {

        final ChimericAlignment.DistancesBetweenAlignmentsOnRefAndOnRead distances = new ChimericAlignment.DistancesBetweenAlignmentsOnRefAndOnRead(simpleChimera);
        final int distBetweenAlignRegionsOnRef = distances.distBetweenAlignRegionsOnRef, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                  distBetweenAlignRegionsOnCtg = distances.distBetweenAlignRegionsOnCtg; // distance-1 between the two regions on contig, denoted as d2 in the comments below
        if ( distBetweenAlignRegionsOnRef > 0 ) {        // Deletion:
            if (distBetweenAlignRegionsOnCtg <= 0) {     // simple deletion when == 0; with homology when < 0
                return TypeInferredFromSimpleChimera.SIMPLE_DEL;
            } else {
                return TypeInferredFromSimpleChimera.RPL;
            }
        } else if (distBetweenAlignRegionsOnRef < 0) {
            if (distBetweenAlignRegionsOnCtg >= 0) { // Tandem repeat expansion:   reference bases [r1e-|d1|+1, r1e] to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1] with optional inserted sequence [c1e+1, c2b-1] in between the two intervals on contig
                return TypeInferredFromSimpleChimera.SMALL_DUP_EXPANSION;
            } else {  // complicated case, see below
                // Deletion:  duplication with repeat number N1 on reference, N2 on contig, such that N1 <= 2*N2 (and N2<N1);
                // Insertion: duplication with repeat number N1 on reference, N2 on contig, such that N2 <= 2*N1 (and N1<N2);
                // in both cases, the equal sign on the right can be taken only when there's pseudo-homology between starting bases of the duplicated sequence and starting bases of the right flanking region
                // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
                return TypeInferredFromSimpleChimera.SMALL_DUP_CPX;
            }
        } else {  // distBetweenAlignRegionsOnRef == 0
            if (distBetweenAlignRegionsOnCtg > 0) { // Insertion: simple insertion, inserted sequence is the sequence [c1e+1, c2b-1] on the contig
                return TypeInferredFromSimpleChimera.SIMPLE_INS;
            } else if (distBetweenAlignRegionsOnCtg < 0) { // Tandem repeat contraction: reference has two copies but one copy was deleted on the contig; duplicated sequence on reference are [r1e-|d2|+1, r1e] and [r2b, r2b+|d2|-1]
                return TypeInferredFromSimpleChimera.DEL_DUP_CONTRACTION;
            } else { // both == 0 => SNP & indel
                throw new GATKException("Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: " + simpleChimera.toString());
            }
        }
    }

    ///////////////

    /**
     * Leaf classes implement their own complication-resolving logic.
     */
    abstract BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                          final byte[] contigSequence,
                                                          final SAMSequenceDictionary referenceDictionary);
    // =================================================================================================================

    /**
     * Intended to be used for simple insertion, simple deletion, and replacement.
     * NOT FOR SMALL DUPLICATIONS!
     *
     * Note that complication resolving logic is to be implemented by leaf classes: simple insertion, simple deletion, and replacement.
     */
    abstract static class SimpleInsertionDeletionBreakpointsInference extends BreakpointsInference {

        protected SimpleInsertionDeletionBreakpointsInference (final ChimericAlignment simpleChimera,
                                                               final byte[] contigSequence,
                                                               final SAMSequenceDictionary referenceDictionary) {

            final BreakpointComplications complications = resolveComplications(simpleChimera, contigSequence, referenceDictionary);

            upstreamBreakpointRefContig
                    = downstreamBreakpointRefContig
                    = simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig();

            final SimpleInterval leftReferenceInterval, rightReferenceInterval;
            if (simpleChimera.isForwardStrandRepresentation) {
                leftReferenceInterval  = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
                rightReferenceInterval = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;
            } else {
                leftReferenceInterval  = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;
                rightReferenceInterval = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
            }
            final int homologyLen = complications.getHomologyForwardStrandRep().length();
            upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
            downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;

            if( complications.insertedSequenceForwardStrandRep.isEmpty() ) { // simple deletion
                altHaplotypeSequence = new byte[0]; // simple deletion has no bases available: (POS, END] is gone!
            } else {
                altHaplotypeSequence = complications.insertedSequenceForwardStrandRep.getBytes();
            }
        }
    }

    ///////////////
    final static class SimpleDeletionBreakpointsInference extends SimpleInsertionDeletionBreakpointsInference {

        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications(simpleChimera, contigSequence);
        }

        SimpleDeletionBreakpointsInference(final ChimericAlignment simpleChimera,
                                           final byte[] contigSequence,
                                           final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);
        }
    }

    /**
     * For class of events where we have assembled across the whole inserted sequence, and with flanking reference sequences.
     */
    final static class InsertionBreakpointsInference extends SimpleInsertionDeletionBreakpointsInference {
        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications(simpleChimera, contigSequence);
        }

        InsertionBreakpointsInference(final ChimericAlignment simpleChimera,
                                      final byte[] contigSequence,
                                      final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);
        }
    }

    final static class ReplacementBreakpointsInference extends SimpleInsertionDeletionBreakpointsInference {
        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return new BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications(simpleChimera, contigSequence);
        }

        ReplacementBreakpointsInference(final ChimericAlignment simpleChimera,
                                        final byte[] contigSequence,
                                        final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);
        }
    }

    ///////////////
    final static class SmallDuplicationBreakpointsInference extends BreakpointsInference {
        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return new BreakpointComplications.SmallDuplicationBreakpointComplications(simpleChimera, contigSequence);
        }

        SmallDuplicationBreakpointsInference(final ChimericAlignment simpleChimera,
                                             final byte[] contigSequence,
                                             final SAMSequenceDictionary referenceDictionary) {

            final BreakpointComplications.SmallDuplicationBreakpointComplications complication =
                    (BreakpointComplications.SmallDuplicationBreakpointComplications) resolveComplications(simpleChimera, contigSequence, referenceDictionary);

            upstreamBreakpointRefContig
                    = downstreamBreakpointRefContig
                    = simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig();

            final int homologyLen = complication.getHomologyForwardStrandRep().length();

            final SimpleInterval leftReferenceInterval, rightReferenceInterval;
            if (simpleChimera.isForwardStrandRepresentation) {
                leftReferenceInterval  = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
                rightReferenceInterval = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;
            } else {
                leftReferenceInterval  = simpleChimera.regionWithHigherCoordOnContig.referenceSpan;
                rightReferenceInterval = simpleChimera.regionWithLowerCoordOnContig.referenceSpan;
            }
            final int expansionUnitDiff = complication.getDupSeqRepeatNumOnCtg() - complication.getDupSeqRepeatNumOnRef();
            if ( expansionUnitDiff > 0 ) {
                upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen
                        - (expansionUnitDiff) * complication.getDupSeqRepeatUnitRefSpan().size();
            } else {
                upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
            }
            downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;

            // TODO: 1/25/18 then extract alt haplotype sequence (to be implemented)
        }
    }

    ///////////////
    abstract static class BNDTypeBreakpointsInference extends BreakpointsInference {
        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return null;
        }
    }

    final static class IntraChrStrandSwitchBreakpointInference extends BNDTypeBreakpointsInference {
        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return new BreakpointComplications.IntraChrStrandSwitchBreakpointComplications(simpleChimera, contigSequence);
        }

        IntraChrStrandSwitchBreakpointInference(final ChimericAlignment simpleChimera,
                                                final byte[] contigSequence,
                                                final SAMSequenceDictionary referenceDictionary) {
            if (simpleChimera.isLikelyInvertedDuplication()){
                final BreakpointComplications.IntraChrStrandSwitchBreakpointComplications complication =
                        (BreakpointComplications.IntraChrStrandSwitchBreakpointComplications) resolveComplications(simpleChimera, contigSequence, referenceDictionary);

                upstreamBreakpointRefPos = complication.getDupSeqRepeatUnitRefSpan().getStart() - 1;
                downstreamBreakpointRefPos = complication.getDupSeqRepeatUnitRefSpan().getEnd();

                altHaplotypeSequence = extractAltHaplotypeForInvDup(simpleChimera, contigSequence);
            }
            else { // simple strand-switch breakpoint with alignments NOT overlapping on the read/assembly contig
                final BreakpointComplications complication = resolveComplications(simpleChimera, contigSequence, referenceDictionary);

                final int homologyLen = complication.getHomologyForwardStrandRep().length();

                final AlignmentInterval one = simpleChimera.regionWithLowerCoordOnContig,
                                        two = simpleChimera.regionWithHigherCoordOnContig;
                final SimpleInterval leftReferenceInterval, rightReferenceInterval;
                if (simpleChimera.isForwardStrandRepresentation) {
                    leftReferenceInterval  = one.referenceSpan;
                    rightReferenceInterval = two.referenceSpan;
                } else {
                    leftReferenceInterval  = two.referenceSpan;
                    rightReferenceInterval = one.referenceSpan;
                }
                if (simpleChimera.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE){
                    upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                    downstreamBreakpointRefPos = rightReferenceInterval.getEnd();
                } else {
                    upstreamBreakpointRefPos = leftReferenceInterval.getStart() - 1;
                    downstreamBreakpointRefPos = rightReferenceInterval.getStart() + homologyLen - 1;
                }
            }

            if( complications.insertedSequenceForwardStrandRep.isEmpty() ) {
                altHaplotypeSequence = new byte[0];
            } else {
                altHaplotypeSequence = complications.insertedSequenceForwardStrandRep.getBytes();
            }
        }

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

    final static class IntraChrRefOrderSwapBreakpointsInference extends BNDTypeBreakpointsInference {
        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return new BreakpointComplications.IntraChrRefOrderSwapBreakpointComplications(simpleChimera, contigSequence);
        }

        IntraChrRefOrderSwapBreakpointsInference(final ChimericAlignment simpleChimera,
                                                 final byte[] contigSequence,
                                                 final SAMSequenceDictionary referenceDictionary) {

            final BreakpointComplications complication = resolveComplications(simpleChimera, contigSequence, referenceDictionary);

            final int homologyLen = complication.getHomologyForwardStrandRep().length();
            final AlignmentInterval one = simpleChimera.regionWithLowerCoordOnContig,
                                    two = simpleChimera.regionWithHigherCoordOnContig;
            final SimpleInterval leftReferenceInterval, rightReferenceInterval;
            if (simpleChimera.isForwardStrandRepresentation) {
                leftReferenceInterval  = one.referenceSpan;
                rightReferenceInterval = two.referenceSpan;
            } else {
                leftReferenceInterval  = two.referenceSpan;
                rightReferenceInterval = one.referenceSpan;
            }
            upstreamBreakpointRefPos = rightReferenceInterval.getStart();
            downstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;

            if( complications.insertedSequenceForwardStrandRep.isEmpty() ) {
                altHaplotypeSequence = new byte[0];
            } else {
                altHaplotypeSequence = complications.insertedSequenceForwardStrandRep.getBytes();
            }
        }
    }

    /**
     * For computing exact and left-adjusted breakpoint locations of inter-chromosome novel adjacency,
     * with or without strand switch.
     */
    final static class InterChromosomeBreakpointsInference extends BreakpointsInference {
        BreakpointComplications resolveComplications(final ChimericAlignment simpleChimera,
                                                     final byte[] contigSequence, SAMSequenceDictionary referenceDictionary) {
            return new BreakpointComplications.InterChromosomeBreakpointComplications(simpleChimera, contigSequence);
        }

        InterChromosomeBreakpointsInference(final ChimericAlignment simpleChimera,
                                            final byte[] contigSequence,
                                            final SAMSequenceDictionary referenceDictionary){

            final BreakpointComplications complication = resolveComplications(simpleChimera, contigSequence, referenceDictionary);

            determineRefContigs(simpleChimera, referenceDictionary);

            extractRefPositions(simpleChimera, complication, referenceDictionary);

            if( complications.insertedSequenceForwardStrandRep.isEmpty() ) {
                altHaplotypeSequence = new byte[0];
            } else {
                altHaplotypeSequence = complications.insertedSequenceForwardStrandRep.getBytes();
            }
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
