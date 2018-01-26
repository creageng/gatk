package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A helper struct for annotating complications that make the locations represented by its associated
 * {@link NovelAdjacencyAndInferredAltHaptype} a little ambiguous.
 *
 * We currently handel five types for simple chimera induced precise variants:
 * <ul>
 *     <li>
 *         deletions
 *     </li>
 *     <li>
 *         small insertions
 *     </li>
 *     <Li>
 *         replacements (i.e. del and ins at the same location)
 *     </Li>
 *     <li>
 *         BND type, which in turn is turned into several subtypes
 *         <ul>
 *             <li>
 *                 intra-chromosome strand-switch BND's, i.e. inversion breakpoint suspects
 *             </li>
 *             <li>
 *                 intra-chromosome no strand-switch but order swap BND's, i.e. large tandem duplication breakpoint suspects
 *             </li>
 *             <li>
 *                 inter-chromosome BND's, with or without strand switch.
 *             </li>
 *         </ul>
 *     </li>
 * </ul>
 */
public abstract class BreakpointComplications {

    /**
     * '+' strand representations of micro-homology, inserted sequence and duplicated sequence on the reference.
     */
    protected String homologyForwardStrandRep = "";
    protected String insertedSequenceForwardStrandRep = "";

    public String getHomologyForwardStrandRep() {
        return homologyForwardStrandRep;
    }
    public String getInsertedSequenceForwardStrandRep() {
        return insertedSequenceForwardStrandRep;
    }

    /**
     * @return Intended for use in debugging and exception message only.
     */
    @Override
    public String toString() {
        return "homology: " + homologyForwardStrandRep + "\tinserted sequence: " + insertedSequenceForwardStrandRep;
    }

    /**
     * Intended to be override by sub classes when more complications are involved.
     */
    public Map<String, Object> toVariantAttributes() {

        final Map<String, Object> attributeMap = new HashMap<>();

        if ( !getInsertedSequenceForwardStrandRep().isEmpty() ) {
            attributeMap.put(GATKSVVCFConstants.INSERTED_SEQUENCE, getInsertedSequenceForwardStrandRep());
        }

        if ( !getHomologyForwardStrandRep().isEmpty() ) {
            attributeMap.put(GATKSVVCFConstants.HOMOLOGY, getHomologyForwardStrandRep());
            attributeMap.put(GATKSVVCFConstants.HOMOLOGY_LENGTH, getHomologyForwardStrandRep().length());
        }

        return attributeMap;
    }

    // =================================================================================================================
    protected BreakpointComplications() {

    }

//    // For test purposes only to cover the deficiency in the default initializations
//    @VisibleForTesting
//    BreakpointComplications(final String homologyForwardStrandRep, final String insertedSequenceForwardStrandRep,
//                            final boolean hasDuplicationAnnotation, final SimpleInterval dupSeqRepeatUnitRefSpan,
//                            final int dupSeqRepeatNumOnRef, final int dupSeqRepeatNumOnCtg,
//                            final List<Strand> dupSeqStrandOnRef, final List<Strand> dupSeqStrandOnCtg,
//                            final List<String> cigarStringsForDupSeqOnCtg, final boolean dupAnnotIsFromOptimization,
//                            final SimpleInterval invertedTransInsertionRefSpan) {
//        this.homologyForwardStrandRep = homologyForwardStrandRep;
//        this.insertedSequenceForwardStrandRep = insertedSequenceForwardStrandRep;
//        this.hasDuplicationAnnotation = hasDuplicationAnnotation;
//        this.dupSeqRepeatUnitRefSpan = dupSeqRepeatUnitRefSpan;
//        this.dupSeqRepeatNumOnRef = dupSeqRepeatNumOnRef;
//        this.dupSeqRepeatNumOnCtg = dupSeqRepeatNumOnCtg;
//        this.dupSeqStrandOnRef = dupSeqStrandOnRef;
//        this.dupSeqStrandOnCtg = dupSeqStrandOnCtg;
//        this.cigarStringsForDupSeqOnCtg = cigarStringsForDupSeqOnCtg;
//        this.dupAnnotIsFromOptimization = dupAnnotIsFromOptimization;
//        this.invertedTransInsertionRefSpan = invertedTransInsertionRefSpan;
//    }

    protected BreakpointComplications(final Kryo kryo, final Input input) {
        homologyForwardStrandRep = input.readString();
        insertedSequenceForwardStrandRep = input.readString();
    }
    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(homologyForwardStrandRep);
        output.writeString(insertedSequenceForwardStrandRep);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BreakpointComplications that = (BreakpointComplications) o;

        if (!homologyForwardStrandRep.equals(that.homologyForwardStrandRep)) return false;
        return insertedSequenceForwardStrandRep.equals(that.insertedSequenceForwardStrandRep);
    }

    @Override
    public int hashCode() {
        int result = homologyForwardStrandRep.hashCode();
        result = 31 * result + insertedSequenceForwardStrandRep.hashCode();
        return result;
    }

    //==================================================================================================================

    /**
     * @return Micro-homology sequence using two alignments of the same contig: as indicated by their overlap on the contig itself.
     *          Empty if they don't overlap on the contig.
     */
    @VisibleForTesting
    static String inferHomology(final AlignmentInterval current, final AlignmentInterval next, final byte[] contigSequence) {

        if (current.endInAssembledContig >= next.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(contigSequence,
                    next.startInAssembledContig-1, current.endInAssembledContig);
            if (current.referenceSpan.getStart() > next.referenceSpan.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            return new String(homologyBytes);
        } else {
            return "";
        }
    }

    /**
     * Note: not suitable for the most complicated case dealt with in {@link BreakpointComplications( ChimericAlignment )}
     * @return Inserted sequence using two alignments of the same contig: as indicated by their separation on the the contig itself.
     */
    @VisibleForTesting
    static String inferInsertedSequence(final AlignmentInterval current, final AlignmentInterval next, final byte[] contigSequence) {

        if (current.endInAssembledContig < next.startInAssembledContig - 1) {
            final byte[] insertedSequenceBytes = Arrays.copyOfRange(contigSequence,
                    current.endInAssembledContig, next.startInAssembledContig - 1);
            if (current.referenceSpan.getStart() > next.referenceSpan.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            return new String(insertedSequenceBytes);
        } else {
            return "";
        }
    }

    //==================================================================================================================

    @DefaultSerializer(SimpleInsDelOrReplacementBreakpointComplications.Serializer.class)
    static final class SimpleInsDelOrReplacementBreakpointComplications extends BreakpointComplications {

        SimpleInsDelOrReplacementBreakpointComplications(final ChimericAlignment simpleChimera, final byte[] contigSeq) {

            final ChimericAlignment.DistancesBetweenAlignmentsOnRefAndOnRead distances = new ChimericAlignment.DistancesBetweenAlignmentsOnRefAndOnRead(simpleChimera);
            final int distBetweenAlignRegionsOnRef = distances.distBetweenAlignRegionsOnRef, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                      distBetweenAlignRegionsOnCtg = distances.distBetweenAlignRegionsOnCtg; // distance-1 between the two regions on contig, denoted as d2 in the comments below

            if ( distBetweenAlignRegionsOnRef > 0 ) { // some bases deleted: here it could be a simple deletion or replacement
                if (distBetweenAlignRegionsOnCtg>=0) {
                    // either: a clean deletion, deleted sequence is [r1e+1, r2b-1] on the reference
                    // or    : deletion with scar, i.e. large non-conserved substitution, reference bases [r1e+1, r2b-1] is substituted with contig bases [c1e+1, c2b-1]
                    insertedSequenceForwardStrandRep = inferInsertedSequence(simpleChimera.regionWithLowerCoordOnContig, simpleChimera.regionWithHigherCoordOnContig, contigSeq);
                } else {
                    // a sequence of bases of length d1+HOM is deleted, and there's homology (which could be dup, but cannot tell): leftFlank+HOM+[r1e+1, r2b-1]+HOM+rightFlank -> leftFlank+HOM+rightFlank
                    homologyForwardStrandRep = inferHomology(simpleChimera.regionWithLowerCoordOnContig, simpleChimera.regionWithHigherCoordOnContig, contigSeq);
                }
            } else if (distBetweenAlignRegionsOnRef == 0 && distBetweenAlignRegionsOnCtg > 0) { // Insertion: simple insertion, inserted sequence is the sequence [c1e+1, c2b-1] on the contig
                insertedSequenceForwardStrandRep = inferInsertedSequence(simpleChimera.regionWithLowerCoordOnContig, simpleChimera.regionWithHigherCoordOnContig, contigSeq);
            } else {
                throw new GATKException("Inferring breakpoint complications with the wrong unit: using simple ins-del unit for simple chimera:\n"
                        + simpleChimera.toString());
            }

            if ( insertedSequenceForwardStrandRep.isEmpty() ){
                    throw new GATKException("An identified breakpoint pair seem to suggest insertion but the inserted sequence is empty: " +
                            simpleChimera.toString());
            }
        }

        protected SimpleInsDelOrReplacementBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleInsDelOrReplacementBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final SimpleInsDelOrReplacementBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public SimpleInsDelOrReplacementBreakpointComplications read(final Kryo kryo, final Input input, final Class<SimpleInsDelOrReplacementBreakpointComplications> klass ) {
                return new SimpleInsDelOrReplacementBreakpointComplications(kryo, input);
            }
        }

        @Override
        public int hashCode() {
            return super.hashCode();
        }

        @Override
        public boolean equals(Object obj) {
            return super.equals(obj);
        }
    }

    @DefaultSerializer(SmallDuplicationBreakpointComplications.Serializer.class)
    static final class SmallDuplicationBreakpointComplications extends BreakpointComplications {
        public static final List<String> DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG = Collections.emptyList();

        private SimpleInterval dupSeqRepeatUnitRefSpan = null;
        private int dupSeqRepeatNumOnRef = 0;
        private int dupSeqRepeatNumOnCtg = 0;
        private List<Strand> dupSeqStrandOnRef = null;
        private List<Strand> dupSeqStrandOnCtg = null;
        private List<String> cigarStringsForDupSeqOnCtg = null;
        private boolean dupAnnotIsFromOptimization = false;

        public SimpleInterval getDupSeqRepeatUnitRefSpan() {
            return dupSeqRepeatUnitRefSpan;
        }
        public int getDupSeqRepeatNumOnRef() {
            return dupSeqRepeatNumOnRef;
        }
        public int getDupSeqRepeatNumOnCtg() {
            return dupSeqRepeatNumOnCtg;
        }
        public List<Strand> getDupSeqOrientationsOnCtg() {
            return dupSeqStrandOnCtg;
        }
        public List<String> getCigarStringsForDupSeqOnCtg() {
            return cigarStringsForDupSeqOnCtg;
        }
        public boolean isDupAnnotIsFromOptimization() {
            return dupAnnotIsFromOptimization;
        }

        @Override
        public Map<String, Object> toVariantAttributes() {
            final Map<String, Object> parentAttributesToBeFilled = super.toVariantAttributes();

            parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, getDupSeqRepeatUnitRefSpan().toString());

            if ( !getCigarStringsForDupSeqOnCtg().isEmpty() ) {
                parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_SEQ_CIGARS,
                        StringUtils.join(getCigarStringsForDupSeqOnCtg(), VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
            }
            parentAttributesToBeFilled.put(GATKSVVCFConstants.DUPLICATION_NUMBERS,
                    new int[]{getDupSeqRepeatNumOnRef(), getDupSeqRepeatNumOnCtg()});
            if ( isDupAnnotIsFromOptimization() ) {
                parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, "");
            }

            if ( getDupSeqOrientationsOnCtg() != null ) {
                parentAttributesToBeFilled.put(GATKSVVCFConstants.DUP_ORIENTATIONS,
                        getDupSeqOrientationsOnCtg().stream().map(Strand::toString).collect(Collectors.joining()));
            }
            return parentAttributesToBeFilled;
        }

        @Override
        public final String toString() {
            String toPrint = super.toString();
            toPrint += String.format("\ttandem duplication repeat unit ref span: %s\t"+
                            "ref repeat num: %d\t"+
                            "ctg repeat num: %d\t"+
                            "dupSeqStrandOnRef: %s\t" +
                            "dupSeqStrandOnCtg: %s\t" +
                            "cigarStringsForDupSeqOnCtg: %s\t"+
                            "tandupAnnotationIsFromSimpleOptimization: %s\t" +
                            "invertedTransInsertionRefSpan: %s",
                    dupSeqRepeatUnitRefSpan == null ? "" : dupSeqRepeatUnitRefSpan,
                    dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg,
                    dupSeqStrandOnRef == null ? "" : dupSeqStrandOnRef.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnRef.size())).toString(),
                    dupSeqStrandOnCtg == null ? "" : dupSeqStrandOnCtg.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnCtg.size())).toString(),
                    cigarStringsForDupSeqOnCtg == null ? "" : cigarStringsForDupSeqOnCtg,
                    isDupAnnotIsFromOptimization() ? "true" : "false");
            return toPrint;
        }


        SmallDuplicationBreakpointComplications(final ChimericAlignment simpleChimera, final byte[] contigSeq) {
            Utils.validateArg(simpleChimera.isNeitherSimpleTranslocationNorIncompletePicture(),
                    "Assumption that the simple chimera is neither incomplete picture nor simple translocation is violated.\n" +
                            simpleChimera.toString());
            Utils.validateArg(simpleChimera.strandSwitch.equals(StrandSwitch.NO_SWITCH),
                    "Assumption that the simple chimera is neither incomplete picture nor simple translocation is violated.\n" +
                            simpleChimera.toString());
            final AlignmentInterval firstContigRegion  = simpleChimera.regionWithLowerCoordOnContig;
            final AlignmentInterval secondContigRegion = simpleChimera.regionWithHigherCoordOnContig;
            final SimpleInterval leftReferenceSpan, rightReferenceSpan;
            if (simpleChimera.isForwardStrandRepresentation) {
                leftReferenceSpan = firstContigRegion.referenceSpan;
                rightReferenceSpan = secondContigRegion.referenceSpan;
            } else {
                leftReferenceSpan = secondContigRegion.referenceSpan;
                rightReferenceSpan = firstContigRegion.referenceSpan;
            }

            final int r1e = leftReferenceSpan.getEnd(),
                    r2b = rightReferenceSpan.getStart(),
                    c1e = firstContigRegion.endInAssembledContig,
                    c2b = secondContigRegion.startInAssembledContig;

            final int distBetweenAlignRegionsOnRef = r2b - r1e - 1, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                      distBetweenAlignRegionsOnCtg = c2b - c1e - 1; // distance-1 between the two regions on contig, denoted as d2 in the comments below

            if (distBetweenAlignRegionsOnRef == 0 && distBetweenAlignRegionsOnCtg < 0) { // Tandem repeat contraction: reference has two copies but one copy was deleted on the contig; duplicated sequence on reference are [r1e-|d2|+1, r1e] and [r2b, r2b+|d2|-1]
                resolveComplicationForSimpleTandupContraction(leftReferenceSpan, firstContigRegion, secondContigRegion, r1e, c1e, c2b, contigSeq);
            } else if (distBetweenAlignRegionsOnRef < 0 && distBetweenAlignRegionsOnCtg >= 0) { // Tandem repeat expansion:   reference bases [r1e-|d1|+1, r1e] to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1] with optional inserted sequence [c1e+1, c2b-1] in between the two intervals on contig
                resolveComplicationForSimpleTandupExpansion(leftReferenceSpan, firstContigRegion, secondContigRegion, r1e, r2b, distBetweenAlignRegionsOnCtg, contigSeq);
            } else if (distBetweenAlignRegionsOnRef < 0 && distBetweenAlignRegionsOnCtg < 0) {  // most complicated case, see below
                // Deletion:  duplication with repeat number N1 on reference, N2 on contig, such that N1 <= 2*N2 (and N2<N1);
                // Insertion: duplication with repeat number N1 on reference, N2 on contig, such that N2 <= 2*N1 (and N1<N2);
                // in both cases, the equal sign on the right can be taken only when there's pseudo-homology between starting bases of the duplicated sequence and starting bases of the right flanking region
                // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
                resolveComplicationForComplexTandup(firstContigRegion, secondContigRegion, r1e, distBetweenAlignRegionsOnRef, distBetweenAlignRegionsOnCtg, contigSeq);
            } else if (distBetweenAlignRegionsOnRef == 0 && distBetweenAlignRegionsOnCtg == 0) {// SNP & indel
                throw new GATKException("Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: " + simpleChimera.toString());
            }

            if ( insertedSequenceForwardStrandRep.isEmpty() ){
                if ( dupSeqRepeatNumOnCtg != dupSeqRepeatNumOnRef && null == dupSeqRepeatUnitRefSpan )
                    throw new GATKException("An identified breakpoint pair seem to suggest insertion but the inserted sequence is empty: " + simpleChimera.toString());
            }
        }

        private void resolveComplicationForSimpleTandupExpansion(final SimpleInterval leftReferenceInterval,
                                                                 final AlignmentInterval firstContigRegion,
                                                                 final AlignmentInterval secondContigRegion,
                                                                 final int r1e, final int r2b,
                                                                 final int distBetweenAlignRegionsOnCtg, final byte[] contigSeq) {
            // note this does not incorporate the duplicated reference sequence
            insertedSequenceForwardStrandRep = distBetweenAlignRegionsOnCtg == 0 ? "" : inferInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
            dupSeqRepeatUnitRefSpan   = new SimpleInterval(leftReferenceInterval.getContig(), r2b, r1e);
            dupSeqRepeatNumOnRef      = 1;
            dupSeqRepeatNumOnCtg      = 2;
            dupSeqStrandOnRef         = Collections.singletonList(Strand.POSITIVE);
            dupSeqStrandOnCtg         = Arrays.asList(Strand.POSITIVE, Strand.POSITIVE);
            cigarStringsForDupSeqOnCtg = new ArrayList<>(2);
            if (firstContigRegion.forwardStrand) {
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(extractCigarForTandupExpansion(firstContigRegion, r1e, r2b)) );
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(extractCigarForTandupExpansion(secondContigRegion, r1e, r2b)) );
            } else {
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandupExpansion(firstContigRegion, r1e, r2b))) );
                cigarStringsForDupSeqOnCtg.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandupExpansion(secondContigRegion, r1e, r2b))) );
            }
        }

        /**
         * Given a {@link AlignmentInterval} from a pair of ARs that forms a {@link ChimericAlignment} signalling a tandem duplication,
         * extract a CIGAR from the {@link AlignmentInterval#cigarAlong5to3DirectionOfContig}
         * that corresponds to the alignment between the suspected repeated sequence on reference between
         * [{@code alignmentIntervalTwoReferenceIntervalSpanBegin}, {@code alignmentIntervalOneReferenceIntervalSpanEnd}],
         * and the sequence in {@link AlignmentInterval#referenceSpan}.
         */
        @VisibleForTesting
        static Cigar extractCigarForTandupExpansion(final AlignmentInterval contigRegion,
                                                    final int alignmentIntervalOneReferenceIntervalSpanEnd,
                                                    final int alignmentIntervalTwoReferenceIntervalSpanBegin) {

            final List<CigarElement> elementList = contigRegion.cigarAlong5to3DirectionOfContig.getCigarElements();
            final List<CigarElement> result = new ArrayList<>(elementList.size());
            final int refStart = contigRegion.referenceSpan.getStart(),
                    refEnd = contigRegion.referenceSpan.getEnd();
            final boolean isForwardStrand = contigRegion.forwardStrand;
            boolean initiatedCollection = false;
            int refPos = isForwardStrand ? refStart : refEnd;
            for(final CigarElement cigarElement : elementList) {
                final CigarOperator operator = cigarElement.getOperator();
                if ( !operator.isClipping() ) {
                    final int opLen = cigarElement.getLength();
                    refPos += operator.consumesReferenceBases() ? (isForwardStrand ? opLen : -opLen) : 0;
                    final int offsetIntoRepeatRegion = isForwardStrand ? refPos - alignmentIntervalTwoReferenceIntervalSpanBegin
                                                                       : alignmentIntervalOneReferenceIntervalSpanEnd - refPos;
                    final int overshootOutOfRepeatRegion = isForwardStrand ? refPos - alignmentIntervalOneReferenceIntervalSpanEnd - 1
                                                                           : alignmentIntervalTwoReferenceIntervalSpanBegin - refPos - 1;

                    if ( offsetIntoRepeatRegion > 0 ) {
                        if ( overshootOutOfRepeatRegion <= 0 ) {
                            result.add( initiatedCollection ? cigarElement : new CigarElement(offsetIntoRepeatRegion, operator));
                            initiatedCollection = true;
                        } else {
                            result.add(new CigarElement(opLen-overshootOutOfRepeatRegion, operator));
                            break;
                        }
                    }
                }
            }

            return new Cigar(result);
        }

        private void resolveComplicationForSimpleTandupContraction(final SimpleInterval leftReferenceInterval,
                                                                   final AlignmentInterval firstContigRegion,
                                                                   final AlignmentInterval secondContigRegion,
                                                                   final int r1e, final int c1e, final int c2b,
                                                                   final byte[] contigSeq) {
            homologyForwardStrandRep = inferHomology(firstContigRegion, secondContigRegion, contigSeq);
            dupSeqRepeatUnitRefSpan  = new SimpleInterval(leftReferenceInterval.getContig(), r1e - ( c1e - c2b ), r1e);
            dupSeqRepeatNumOnRef     = 2;
            dupSeqRepeatNumOnCtg     = 1;
            dupSeqStrandOnRef        = Arrays.asList(Strand.POSITIVE, Strand.POSITIVE);
            dupSeqStrandOnCtg        = Collections.singletonList(Strand.POSITIVE);
            cigarStringsForDupSeqOnCtg = DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG;
        }

        private void resolveComplicationForComplexTandup(final AlignmentInterval firstContigRegion,
                                                         final AlignmentInterval secondContigRegion,
                                                         final int r1e, final int distBetweenAlignRegionsOnRef,
                                                         final int distBetweenAlignRegionsOnCtg, final byte[] contigSeq) {

            final TandemRepeatStructure duplicationComplication =
                    new TandemRepeatStructure(distBetweenAlignRegionsOnRef, distBetweenAlignRegionsOnCtg);

            final boolean isExpansion     = distBetweenAlignRegionsOnRef<distBetweenAlignRegionsOnCtg;

            final int repeatUnitSpanStart = r1e - duplicationComplication.pseudoHomologyLen
                    - duplicationComplication.repeatedSeqLen * duplicationComplication.lowerRepeatNumberEstimate
                    + 1;
            final int repeatUnitSpanEnd   = repeatUnitSpanStart + duplicationComplication.repeatedSeqLen - 1;
            homologyForwardStrandRep      = inferHomology(firstContigRegion, secondContigRegion, contigSeq);
            cigarStringsForDupSeqOnCtg    = DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG;
            dupSeqRepeatUnitRefSpan       = new SimpleInterval(firstContigRegion.referenceSpan.getContig(), repeatUnitSpanStart, repeatUnitSpanEnd);
            dupSeqRepeatNumOnRef          = isExpansion ? duplicationComplication.lowerRepeatNumberEstimate
                    : duplicationComplication.higherRepeatNumberEstimate;
            dupSeqRepeatNumOnCtg          = isExpansion ? duplicationComplication.higherRepeatNumberEstimate
                    : duplicationComplication.lowerRepeatNumberEstimate;
            dupSeqStrandOnRef             = new ArrayList<>(Collections.nCopies(dupSeqRepeatNumOnRef, Strand.POSITIVE));
            dupSeqStrandOnCtg             = new ArrayList<>(Collections.nCopies(dupSeqRepeatNumOnCtg, Strand.POSITIVE));
            dupAnnotIsFromOptimization    = true;
        }

        // TODO: 03/03/17 this complicated tandem duplication annotation is not exactly reproducible in the following sense:
        //          1) depending on what the assembler might produce, e.g. different runs producing slightly different sequences
        //          hence affecting alignment,
        //          2) the assembler might decide to output RC sequences between runs hence the mapping would be to '+' or '-' strand
        //       these randomness may give slightly different results by this treatment
        /**
         * This auxiliary structure, when constructed given overlaps of two corresponding regions on reference and contig sequences,
         * attempts to find--naively and slowly--the repeat numbers on the reference and on the contig of tandem repeats,
         * as well as the pseudo-homology between the duplicated sequence and the right flanking region.
         *
         * An example might help:
         * an assembled contig that's actually a repeat expansion from 1 repeat to 2 repeats with pseudo-homology:
         * TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
         * can be aligned to chr18,
         * the 1st alignment chr18:312579-718, 140M135S, which can be broken into the following part
         * 31:  TGCCAGGTTACATGGCAAAGAGGGTAGATAT
         * 109: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAA
         * 135: GAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
         * And the arithmetic to get the cigar operation length works this way:
         * 31 + 109 = 140
         * 109 = 96 + 13
         * where 31 is the left flanking region before the repeated unit, which itself is 96 bases long (see below),
         * the number 13 is the length of the pseudo-homology between the starting bases of the repeated sequence and the right flanking region
         * a clearer picture emerges when we look at the 2nd alignment
         * chr18:312610-757, 127S148M, which can be broken into
         * 31: TGCCAGGTTACATGGCAAAGAGGGTAGATAT
         * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA
         * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA
         * 13: GGGCAGCTGTGGA
         * 39: TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
         * And the arithmetic works this way:
         * 31 + 96 = 127
         * 96 + 13 + 39 = 148
         */
        private static final class TandemRepeatStructure {

            /**
             * In {@link TandemRepeatStructure} where the naive attempt to resolve number of tandem repeats
             * on the reference and sample is done, we assume the lower number of repeats is no higher than this number.
             */
            private static final int MAX_LOWER_CN = 10;

            final int lowerRepeatNumberEstimate;
            final int higherRepeatNumberEstimate;
            final int repeatedSeqLen;
            final int pseudoHomologyLen;


            @VisibleForTesting
            TandemRepeatStructure(final int distBetweenAlignRegionsOnRef, final int distBetweenAlignRegionsOnCtg) {
                // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
                final boolean isExpansion = distBetweenAlignRegionsOnRef < distBetweenAlignRegionsOnCtg;
                final int overlapOnLowerCNSequence, overlapOnHigherCNSequence;
                if (isExpansion) {
                    overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
                    overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
                } else {     // d1 is lower absolute value -> reference has higher copy number of the duplication, i.e. Deletion
                    overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
                    overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
                }

                int higherCnEst=0, lowerCnEst=0, unitLen=0, pseudoHomLen=0;
                double err = Double.MAX_VALUE;
                for(int cn2 = 1; cn2< MAX_LOWER_CN; ++cn2) {
                    for(int cn1 = cn2 + 1; cn1 <= 2 * cn2; ++cn1) {
                        final int dupLenUpperBound = (cn1 == 2 * cn2) ? overlapOnLowerCNSequence : overlapOnHigherCNSequence;
                        for (int l = 2; l <= dupLenUpperBound; ++l) {
                            for (int lambda = 0; lambda < l; ++lambda) {
                                final int d1 = (2*cn2 - cn1)*l + lambda;
                                final int d2 = cn2*l + lambda;
                                final double newErr = Math.abs(overlapOnHigherCNSequence-d1) + Math.abs(overlapOnLowerCNSequence-d2);
                                if (newErr < err) {
                                    err = newErr;
                                    higherCnEst = cn1; lowerCnEst = cn2;
                                    unitLen= l; pseudoHomLen = lambda;
                                }
                                if (err < 1){
                                    lowerRepeatNumberEstimate = lowerCnEst;
                                    higherRepeatNumberEstimate = higherCnEst;
                                    repeatedSeqLen = unitLen;
                                    pseudoHomologyLen = pseudoHomLen;
                                    return;
                                }
                            }
                        }
                    }
                }

                lowerRepeatNumberEstimate = lowerCnEst;
                higherRepeatNumberEstimate = higherCnEst;
                repeatedSeqLen = unitLen;
                pseudoHomologyLen = pseudoHomLen;
            }
        }

        protected SmallDuplicationBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
            final String ctg = input.readString();
            final int start = input.readInt();
            final int end = input.readInt();
            dupSeqRepeatUnitRefSpan = new SimpleInterval(ctg, start, end);
            dupSeqRepeatNumOnRef = input.readInt();
            dupSeqRepeatNumOnCtg = input.readInt();
            dupSeqStrandOnRef = new ArrayList<>(dupSeqRepeatNumOnRef);
            for (int i=0; i<dupSeqRepeatNumOnRef; ++i) {
                dupSeqStrandOnRef.add(Strand.values()[input.readInt()]);
            }
            dupSeqStrandOnCtg = new ArrayList<>(dupSeqRepeatNumOnCtg);
            for (int i=0; i<dupSeqRepeatNumOnCtg; ++i) {
                dupSeqStrandOnCtg.add(Strand.values()[input.readInt()]);
            }
            final int cigarCounts = input.readInt();
            cigarStringsForDupSeqOnCtg = new ArrayList<>(cigarCounts);
            for(int i = 0; i < cigarCounts; ++i) {
                cigarStringsForDupSeqOnCtg.add(input.readString());
            }
            dupAnnotIsFromOptimization = input.readBoolean();
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
            output.writeString(dupSeqRepeatUnitRefSpan.getContig());
            output.writeInt(dupSeqRepeatUnitRefSpan.getStart());
            output.writeInt(dupSeqRepeatUnitRefSpan.getEnd());
            output.writeInt(dupSeqRepeatNumOnRef);
            output.writeInt(dupSeqRepeatNumOnCtg);
            dupSeqStrandOnRef.forEach(s -> output.writeInt(s.ordinal()));
            dupSeqStrandOnCtg.forEach(s -> output.writeInt(s.ordinal()));
            output.writeInt(cigarStringsForDupSeqOnCtg.size());
            cigarStringsForDupSeqOnCtg.forEach(output::writeString);
            output.writeBoolean(dupAnnotIsFromOptimization);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SmallDuplicationBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final SmallDuplicationBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public SmallDuplicationBreakpointComplications read(final Kryo kryo, final Input input, final Class<SmallDuplicationBreakpointComplications> klass ) {
                return new SmallDuplicationBreakpointComplications(kryo, input);
            }
        }
    }

    abstract static class BNDTypeBreakpointComplications extends BreakpointComplications {
        protected BNDTypeBreakpointComplications(final ChimericAlignment simpleChimera, final byte[] contigSeq) {
            homologyForwardStrandRep = inferHomology(simpleChimera.regionWithLowerCoordOnContig,
                                                    simpleChimera.regionWithHigherCoordOnContig, contigSeq);
            insertedSequenceForwardStrandRep = inferInsertedSequence(simpleChimera.regionWithLowerCoordOnContig,
                                                                    simpleChimera.regionWithHigherCoordOnContig, contigSeq);
        }

        protected BNDTypeBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }
    }

    /**
     * For novel adjacency between reference locations that are on the same chromosome, and with a strand switch
     */
    @DefaultSerializer(IntraChrStrandSwitchBreakpointComplications.Serializer.class)
    static final class IntraChrStrandSwitchBreakpointComplications extends BNDTypeBreakpointComplications {
        public static final List<String> DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG = Collections.emptyList();

        private SimpleInterval dupSeqRepeatUnitRefSpan = null;
        private int dupSeqRepeatNumOnRef = 0;
        private int dupSeqRepeatNumOnCtg = 0;
        private List<Strand> dupSeqStrandOnRef = null;
        private List<Strand> dupSeqStrandOnCtg = null;
        private List<String> cigarStringsForDupSeqOnCtg = null;
        private boolean dupAnnotIsFromOptimization = false;

        public SimpleInterval getDupSeqRepeatUnitRefSpan() {
            return dupSeqRepeatUnitRefSpan;
        }
        public int getDupSeqRepeatNumOnRef() {
            return dupSeqRepeatNumOnRef;
        }
        public int getDupSeqRepeatNumOnCtg() {
            return dupSeqRepeatNumOnCtg;
        }
        public List<Strand> getDupSeqOrientationsOnCtg() {
            return dupSeqStrandOnCtg;
        }
        public List<String> getCigarStringsForDupSeqOnCtg() {
            return cigarStringsForDupSeqOnCtg;
        }
        public boolean isDupAnnotIsFromOptimization() {
            return dupAnnotIsFromOptimization;
        }

        static final List<Strand> DEFAULT_INV_DUP_REF_ORIENTATION = Collections.singletonList(Strand.POSITIVE);
        static final List<Strand> DEFAULT_INV_DUP_CTG_ORIENTATIONS_FR = Arrays.asList(Strand.POSITIVE, Strand.NEGATIVE);
        static final List<Strand> DEFAULT_INV_DUP_CTG_ORIENTATIONS_RF = Arrays.asList(Strand.NEGATIVE, Strand.POSITIVE);
        private SimpleInterval invertedTransInsertionRefSpan = null; // TODO: 10/2/17 see ticket 3647
        public SimpleInterval getInvertedTransInsertionRefSpan() {
            return invertedTransInsertionRefSpan;
        }

        @Override
        public final String toString() {
            final String toPrint = super.toString();
            toPrint += String.format("\ttandem duplication repeat unit ref span: %s\t"+
                            "ref repeat num: %d\t"+
                            "ctg repeat num: %d\t"+
                            "dupSeqStrandOnRef: %s\t" +
                            "dupSeqStrandOnCtg: %s\t" +
                            "cigarStringsForDupSeqOnCtg: %s\t"+
                            "tandupAnnotationIsFromSimpleOptimization: %s\t" +
                            "invertedTransInsertionRefSpan: %s",
                    dupSeqRepeatUnitRefSpan == null ? "" : dupSeqRepeatUnitRefSpan,
                    dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg,
                    dupSeqStrandOnRef == null ? "" : dupSeqStrandOnRef.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnRef.size())).toString(),
                    dupSeqStrandOnCtg == null ? "" : dupSeqStrandOnCtg.stream().map(Strand::toString).collect(SVUtils.arrayListCollector(dupSeqStrandOnCtg.size())).toString(),
                    cigarStringsForDupSeqOnCtg == null ? "" : cigarStringsForDupSeqOnCtg,
                    isDupAnnotIsFromOptimization() ? "true" : "false",
                    invertedTransInsertionRefSpan == null ? "" : invertedTransInsertionRefSpan);
            return toPrint;
        }

        IntraChrStrandSwitchBreakpointComplications(final ChimericAlignment simpleChimera, final byte[] contigSeq) {
            super(simpleChimera, contigSeq);
            if ( simpleChimera.isLikelyInvertedDuplication() ) {
                resolveComplicationForInvDup(simpleChimera, contigSeq);
            } else {
                resolveComplicationForSimpleStrandSwitch(simpleChimera, contigSeq);
            }
        }

        void resolveComplicationForSimpleStrandSwitch(final ChimericAlignment chimericAlignment, final byte[] contigSeq) {

            final AlignmentInterval firstAlignmentInterval  = chimericAlignment.regionWithLowerCoordOnContig;
            final AlignmentInterval secondAlignmentInterval = chimericAlignment.regionWithHigherCoordOnContig;

            homologyForwardStrandRep = inferHomology(firstAlignmentInterval, secondAlignmentInterval, contigSeq);
            insertedSequenceForwardStrandRep = inferInsertedSequence(firstAlignmentInterval, secondAlignmentInterval, contigSeq);
            dupSeqRepeatUnitRefSpan = null;
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
            dupSeqStrandOnRef = dupSeqStrandOnCtg = null;
            cigarStringsForDupSeqOnCtg = null;
            dupAnnotIsFromOptimization = false;
            hasDuplicationAnnotation = false;
        }

        /**
         * Initialize the fields in this object, assuming the input chimeric alignment is induced by two alignments with
         * "significant" (see {@link ChimericAlignment#isLikelyInvertedDuplication()})
         * overlap on their reference spans.
         */
        private void resolveComplicationForInvDup(final ChimericAlignment chimericAlignment, final byte[] contigSeq) {

            final AlignmentInterval firstAlignmentInterval  = chimericAlignment.regionWithLowerCoordOnContig;
            final AlignmentInterval secondAlignmentInterval = chimericAlignment.regionWithHigherCoordOnContig;

            // TODO: 8/8/17 this might be wrong regarding how strand is involved, fix it
            insertedSequenceForwardStrandRep = inferInsertedSequence(firstAlignmentInterval, secondAlignmentInterval, contigSeq);
            hasDuplicationAnnotation = true;

            dupSeqRepeatNumOnRef = 1;
            dupSeqRepeatNumOnCtg = 2;
            dupSeqStrandOnRef = DEFAULT_INV_DUP_REF_ORIENTATION;

            // jump start and jump landing locations
            final int jumpStartRefLoc = firstAlignmentInterval.forwardStrand ? firstAlignmentInterval.referenceSpan.getEnd()
                    : firstAlignmentInterval.referenceSpan.getStart();
            final int jumpLandingRefLoc = secondAlignmentInterval.forwardStrand ? secondAlignmentInterval.referenceSpan.getStart()
                    : secondAlignmentInterval.referenceSpan.getEnd();

            if (firstAlignmentInterval.forwardStrand) {
                final int alpha = firstAlignmentInterval.referenceSpan.getStart(),
                        omega = secondAlignmentInterval.referenceSpan.getStart();
                dupSeqRepeatUnitRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                        Math.max(alpha, omega), Math.min(jumpStartRefLoc, jumpLandingRefLoc));
                if ( (alpha <= omega && jumpStartRefLoc < jumpLandingRefLoc) || (alpha > omega && jumpLandingRefLoc < jumpStartRefLoc) ) {
                    invertedTransInsertionRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                            Math.min(jumpStartRefLoc, jumpLandingRefLoc) + 1, Math.max(jumpStartRefLoc, jumpLandingRefLoc));
                }
                dupSeqStrandOnCtg = DEFAULT_INV_DUP_CTG_ORIENTATIONS_FR;
            } else {
                final int alpha = firstAlignmentInterval.referenceSpan.getEnd(),
                        omega = secondAlignmentInterval.referenceSpan.getEnd();
                dupSeqRepeatUnitRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                        Math.max(jumpStartRefLoc, jumpLandingRefLoc), Math.min(alpha, omega));
                if ( (alpha >= omega && jumpLandingRefLoc < jumpStartRefLoc) || (alpha < omega && jumpStartRefLoc < jumpLandingRefLoc) ) {
                    invertedTransInsertionRefSpan = new SimpleInterval(firstAlignmentInterval.referenceSpan.getContig(),
                            Math.min(jumpStartRefLoc, jumpLandingRefLoc) + 1, Math.max(jumpStartRefLoc, jumpLandingRefLoc));
                }
                dupSeqStrandOnCtg = DEFAULT_INV_DUP_CTG_ORIENTATIONS_RF;
            }
            cigarStringsForDupSeqOnCtg = DEFAULT_CIGAR_STRINGS_FOR_DUP_SEQ_ON_CTG; // not computing cigars because alt haplotypes will be extracted

            dupAnnotIsFromOptimization = false;
        }

        protected IntraChrStrandSwitchBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);

            if (input.readBoolean()) {
                final String chr = input.readString();
                final int start = input.readInt();
                final int end = input.readInt();
                invertedTransInsertionRefSpan = new SimpleInterval(chr, start, end);
            }
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
            output.writeString(dupSeqRepeatUnitRefSpan.getContig());
            output.writeInt(dupSeqRepeatUnitRefSpan.getStart());
            output.writeInt(dupSeqRepeatUnitRefSpan.getEnd());
            output.writeInt(dupSeqRepeatNumOnRef);
            output.writeInt(dupSeqRepeatNumOnCtg);
            dupSeqStrandOnRef.forEach(s -> output.writeInt(s.ordinal()));
            dupSeqStrandOnCtg.forEach(s -> output.writeInt(s.ordinal()));
            output.writeInt(cigarStringsForDupSeqOnCtg.size());
            cigarStringsForDupSeqOnCtg.forEach(output::writeString);
            output.writeBoolean(dupAnnotIsFromOptimization);

            output.writeBoolean(invertedTransInsertionRefSpan != null);
            if (invertedTransInsertionRefSpan != null) {
                output.writeString(invertedTransInsertionRefSpan.getContig());
                output.writeInt(invertedTransInsertionRefSpan.getStart());
                output.writeInt(invertedTransInsertionRefSpan.getEnd());
            }
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<IntraChrStrandSwitchBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final IntraChrStrandSwitchBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public IntraChrStrandSwitchBreakpointComplications read(final Kryo kryo, final Input input, final Class<IntraChrStrandSwitchBreakpointComplications> klass ) {
                return new IntraChrStrandSwitchBreakpointComplications(kryo, input);
            }
        }
    }

    /**
     * For novel adjacency between reference locations that are on the same chromosome,
     * WITHOUT strand switch but with order swap,
     * i.e. a base with higher coordinate on ref has lower coordinate on sample.
     */
    @DefaultSerializer(IntraChrRefOrderSwapBreakpointComplications.Serializer.class)
    static final class IntraChrRefOrderSwapBreakpointComplications extends BNDTypeBreakpointComplications {

        IntraChrRefOrderSwapBreakpointComplications(final ChimericAlignment simpleChimera, final byte[] contigSeq) {
            super(simpleChimera, contigSeq);
        }

        private IntraChrRefOrderSwapBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);

        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<IntraChrRefOrderSwapBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final IntraChrRefOrderSwapBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public IntraChrRefOrderSwapBreakpointComplications read(final Kryo kryo, final Input input, final Class<IntraChrRefOrderSwapBreakpointComplications> klass ) {
                return new IntraChrRefOrderSwapBreakpointComplications(kryo, input);
            }
        }
    }

    @DefaultSerializer(InterChromosomeBreakpointComplications.Serializer.class)
    static final class InterChromosomeBreakpointComplications extends BNDTypeBreakpointComplications  {
        InterChromosomeBreakpointComplications(final ChimericAlignment simpleChimera, final byte[] contigSeq) {
            super(simpleChimera, contigSeq);
        }

        protected InterChromosomeBreakpointComplications(final Kryo kryo, final Input input) {
            super(kryo, input);
        }

        protected void serialize(final Kryo kryo, final Output output) {
            super.serialize(kryo, output);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<InterChromosomeBreakpointComplications> {
            @Override
            public void write(final Kryo kryo, final Output output, final InterChromosomeBreakpointComplications breakpointComplications) {
                breakpointComplications.serialize(kryo, output);
            }

            @Override
            public InterChromosomeBreakpointComplications read(final Kryo kryo, final Input input, final Class<InterChromosomeBreakpointComplications> klass ) {
                return new InterChromosomeBreakpointComplications(kryo, input);
            }
        }
    }
}
