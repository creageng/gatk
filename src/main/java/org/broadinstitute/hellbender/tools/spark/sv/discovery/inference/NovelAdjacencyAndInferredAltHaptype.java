package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.Arrays;

/**
 * This class represents a pair of inferred genomic locations on the reference whose novel adjacency is generated
 * due to a simple SV event (in other words, a simple rearrangement between two genomic locations)
 * that is suggested by the input {@link AlignedContig},
 * and complications in pinning down the locations to exact base pair resolution.
 */
@DefaultSerializer(NovelAdjacencyAndInferredAltHaptype.Serializer.class)
public class NovelAdjacencyAndInferredAltHaptype {

    public final SimpleInterval leftJustifiedLeftRefLoc;
    public final SimpleInterval leftJustifiedRightRefLoc;

    public final StrandSwitch strandSwitch;
    public final BreakpointComplications complication;

    public final byte[] altHaplotypeSequence;

    public NovelAdjacencyAndInferredAltHaptype(final ChimericAlignment chimericAlignment, final byte[] contigSequence,
                                               final SAMSequenceDictionary referenceDictionary) {

        strandSwitch = chimericAlignment.strandSwitch;

        try {
            complication = new BreakpointComplications(chimericAlignment, contigSequence);

            final Tuple2<SimpleInterval, SimpleInterval> leftJustifiedBreakpoints =
                    leftJustifyBreakpoints(chimericAlignment, complication, referenceDictionary);
            leftJustifiedLeftRefLoc = leftJustifiedBreakpoints._1();
            leftJustifiedRightRefLoc = leftJustifiedBreakpoints._2();

            altHaplotypeSequence =
                    SimpleNovelAdjacencyInterpreter.BreakpointsInference.getInferenceClass(chimericAlignment, complication, referenceDictionary)
                            .getInferredAltHaplotypeSequence();
        } catch (final IllegalArgumentException iaex) { // catching IAEX specifically because it is the most likely exception thrown if there's bug, this helps quickly debugging what the problem is
            throw new GATKException("Erred when inferring breakpoint location and event type from chimeric alignment:\n" +
                    chimericAlignment.toString(), iaex);
        }
    }

    protected NovelAdjacencyAndInferredAltHaptype(final Kryo kryo, final Input input) {
        final String contig1 = input.readString();
        final int start1 = input.readInt();
        final int end1 = input.readInt();
        this.leftJustifiedLeftRefLoc = new SimpleInterval(contig1, start1, end1);
        final String contig2 = input.readString();
        final int start2 = input.readInt();
        final int end2 = input.readInt();
        this.leftJustifiedRightRefLoc = new SimpleInterval(contig2, start2, end2);

        this.strandSwitch = StrandSwitch.values()[input.readInt()];
        this.complication = kryo.readObject(input, BreakpointComplications.class);

        final boolean altSeqIsNull = input.readBoolean();
        if (altSeqIsNull) {
            altHaplotypeSequence = null;
        } else {
            final int arraySize = input.readInt();
            altHaplotypeSequence = new byte[arraySize];
            for (int i = 0 ; i < arraySize; ++i) {
                altHaplotypeSequence[i] = input.readByte();
            }
        }
    }

    /**
     * Returns the reference coordinates of the upstream and downstream breakpoints implied by input {@code chimericAlignment}.
     * If there is homologous sequence represented in the alignments, it will be assigned to the side of the breakpoint
     * with higher reference coordinates, i.e. we follow left alignment convention here.
     */
    @VisibleForTesting
    static Tuple2<SimpleInterval, SimpleInterval> leftJustifyBreakpoints(final ChimericAlignment chimericAlignment,
                                                                         final BreakpointComplications complication,
                                                                         final SAMSequenceDictionary referenceDictionary) {

        final Tuple2<SimpleInterval, SimpleInterval> leftAdjustedBreakpoints =
                SimpleNovelAdjacencyInterpreter.BreakpointsInference.getInferenceClass(chimericAlignment, complication, referenceDictionary)
                        .getLeftJustifiedBreakpoints();

        validateInferredLocations(leftAdjustedBreakpoints._1, leftAdjustedBreakpoints._2, referenceDictionary, chimericAlignment, complication);

        return leftAdjustedBreakpoints;
    }

    private static void validateInferredLocations(final SimpleInterval leftBreakpoint, final SimpleInterval rightBreakpoint,
                                                  final SAMSequenceDictionary referenceSequenceDictionary,
                                                  final ChimericAlignment ca, final BreakpointComplications complication) {

        Utils.validate(IntervalUtils.isBefore(leftBreakpoint, rightBreakpoint, referenceSequenceDictionary) ||
                        leftBreakpoint.equals(rightBreakpoint),
                "Inferred novel adjacency reference locations have left location after right location : left " +
                        leftBreakpoint + "\tright " + rightBreakpoint + "\t"  + ca.toString() + "\n" + complication.toString());

        Utils.validate(leftBreakpoint.getEnd() <= referenceSequenceDictionary.getSequence(leftBreakpoint.getContig()).getSequenceLength(),
                "Inferred breakpoint beyond reference sequence length: inferred location: " + leftBreakpoint +
                        "\tref contig length: " + referenceSequenceDictionary.getSequence(leftBreakpoint.getContig()).getSequenceLength() + "\n"
                        + ca.toString() + "\n" + complication.toString());
        Utils.validate(rightBreakpoint.getEnd() <= referenceSequenceDictionary.getSequence(rightBreakpoint.getContig()).getSequenceLength(),
                "Inferred breakpoint beyond reference sequence length: inferred location: " + rightBreakpoint +
                        "\tref contig length: " + referenceSequenceDictionary.getSequence(rightBreakpoint.getContig()).getSequenceLength() + "\n"
                        + ca.toString() + "\n" + complication.toString());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        NovelAdjacencyAndInferredAltHaptype that = (NovelAdjacencyAndInferredAltHaptype) o;

        if (!leftJustifiedLeftRefLoc.equals(that.leftJustifiedLeftRefLoc)) return false;
        if (!leftJustifiedRightRefLoc.equals(that.leftJustifiedRightRefLoc)) return false;
        if (strandSwitch != that.strandSwitch) return false;
        if (!complication.equals(that.complication)) return false;
        return Arrays.equals(altHaplotypeSequence, that.altHaplotypeSequence);
    }

    @Override
    public int hashCode() {
        int result = leftJustifiedLeftRefLoc.hashCode();
        result = 31 * result + leftJustifiedRightRefLoc.hashCode();
        result = 31 * result + strandSwitch.hashCode();
        result = 31 * result + complication.hashCode();
        result = 31 * result + Arrays.hashCode(altHaplotypeSequence);
        return result;
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftJustifiedLeftRefLoc.getContig());
        output.writeInt(leftJustifiedLeftRefLoc.getStart());
        output.writeInt(leftJustifiedLeftRefLoc.getEnd());
        output.writeString(leftJustifiedRightRefLoc.getContig());
        output.writeInt(leftJustifiedRightRefLoc.getStart());
        output.writeInt(leftJustifiedRightRefLoc.getEnd());
        output.writeInt(strandSwitch.ordinal());
        kryo.writeObject(output, complication);

        if (altHaplotypeSequence==null) {
            output.writeBoolean(true);
        } else {
            output.writeBoolean(false);
            output.writeInt(altHaplotypeSequence.length);
            for (final byte b : altHaplotypeSequence) {
                output.writeByte(b);
            }
        }
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NovelAdjacencyAndInferredAltHaptype> {
        @Override
        public void write(final Kryo kryo, final Output output, final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype) {
            novelAdjacencyAndInferredAltHaptype.serialize(kryo, output);
        }

        @Override
        public NovelAdjacencyAndInferredAltHaptype read(final Kryo kryo, final Input input, final Class<NovelAdjacencyAndInferredAltHaptype> klass ) {
            return new NovelAdjacencyAndInferredAltHaptype(kryo, input);
        }
    }

    /**
     * @return Intended for use in debugging and exception message only.
     */
    @Override
    public String toString() {
        return String.format("%s\t%s\t%s\t%s", leftJustifiedLeftRefLoc.toString(), leftJustifiedRightRefLoc.toString(),
                strandSwitch.name(), complication.toString());
    }
}
