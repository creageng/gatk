package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndInferredAltHaptype;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.Collections;
import java.util.Map;

public abstract class SimpleSVType extends SvType {
    public static String createBracketedSymbAlleleString(final String vcfHeaderDefinedSymbAltAllele) {
        return "<" + vcfHeaderDefinedSymbAltAllele + ">";
    }

    protected SimpleSVType(final String id, final Allele altAllele, final int len, final Map<String, String> typeSpecificExtraAttributes) {
        super(id, altAllele, len, typeSpecificExtraAttributes);
    }

    public enum TYPES {
        INV, DEL, INS, DUP, DUP_INV;
    }

    public static final class Inversion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.INV.name();
        }

        public Inversion(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            super(getIDString(novelAdjacencyAndInferredAltHaptype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INV)),
                    novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart() - novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd(),
                    Collections.singletonMap((novelAdjacencyAndInferredAltHaptype.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33, ""));
        }

        private static String getIDString(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            final String contig = novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getContig();
            final int start = novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd();
            final int end = novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart();
            final StrandSwitch strandSwitch = novelAdjacencyAndInferredAltHaptype.strandSwitch;

            return (strandSwitch.equals(StrandSwitch.FORWARD_TO_REVERSE) ? GATKSVVCFConstants.INV55 : GATKSVVCFConstants.INV33) + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    contig + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + start + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR + end;
        }
    }

    public static final class Deletion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DEL.name();
        }

        @SuppressWarnings("unchecked")
        public Deletion(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            super(getIDString(novelAdjacencyAndInferredAltHaptype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    -(novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart() - novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd()),
                    novelAdjacencyAndInferredAltHaptype.complication.hasDuplicationAnnotation() ? Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING, "") : Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {

            return  ((novelAdjacencyAndInferredAltHaptype.complication.hasDuplicationAnnotation()) ? GATKSVVCFConstants.DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING : TYPES.DEL.name())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart();
        }
    }

    public static final class Insertion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.INS.name();
        }

        @SuppressWarnings("unchecked")
        public Insertion(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            super(getIDString(novelAdjacencyAndInferredAltHaptype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INS)),
                    novelAdjacencyAndInferredAltHaptype.complication.getInsertedSequenceForwardStrandRep().length(),
                    Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {

            return TYPES.INS.name() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart();
        }
    }

    public static final class DuplicationTandem extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DUP.name();
        }

        public DuplicationTandem(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            super(getIDString(novelAdjacencyAndInferredAltHaptype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP)),
                    novelAdjacencyAndInferredAltHaptype.complication.getInsertedSequenceForwardStrandRep().length()
                            + (novelAdjacencyAndInferredAltHaptype.complication.getDupSeqRepeatNumOnCtg() - novelAdjacencyAndInferredAltHaptype.complication.getDupSeqRepeatNumOnRef())* novelAdjacencyAndInferredAltHaptype.complication.getDupSeqRepeatUnitRefSpan().size(),
                    Collections.singletonMap(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, ""));
        }

        private static String getIDString(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {

            return GATKSVVCFConstants.DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart();
        }
    }

    public static final class ImpreciseDeletion extends SimpleSVType {

        @Override
        public String toString() {
            return TYPES.DEL.name();
        }

        @SuppressWarnings("unchecked")
        public ImpreciseDeletion(final EvidenceTargetLink evidenceTargetLink, final ReadMetadata metadata) {

            super(getIDString(evidenceTargetLink, metadata),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL)),
                    (evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().midpoint() -
                            evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval().midpoint()),
                    Collections.EMPTY_MAP);
        }

        private static String getIDString(final EvidenceTargetLink evidenceTargetLink, final ReadMetadata metadata) {

            return TYPES.DEL.name()
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + GATKSVVCFConstants.IMPRECISE + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + metadata.getContigName(evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().getContig())
                    + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getLeft().getInterval().getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval().getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + evidenceTargetLink.getPairedStrandedIntervals().getRight().getInterval().getEnd();
        }
    }

    public static final class DuplicationInverted extends SimpleSVType {

        @Override
        public String toString() {
            return "DUP:INV";
        }

        @SuppressWarnings("unchecked")
        public DuplicationInverted(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            super(getIDString(novelAdjacencyAndInferredAltHaptype),
                    Allele.create(createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_INVDUP)),
                    novelAdjacencyAndInferredAltHaptype.complication.getDupSeqRepeatUnitRefSpan().size(),
                    Collections.EMPTY_MAP);
        }

        private static String getIDString(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            return GATKSVVCFConstants.DUP_INV_INTERNAL_ID_START_STRING + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                    + novelAdjacencyAndInferredAltHaptype.leftJustifiedRightRefLoc.getStart();
        }
    }

}
