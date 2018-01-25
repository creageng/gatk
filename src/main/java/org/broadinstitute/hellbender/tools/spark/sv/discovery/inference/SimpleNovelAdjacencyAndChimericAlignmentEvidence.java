package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

@DefaultSerializer(SimpleNovelAdjacencyAndChimericAlignmentEvidence.Serializer.class)
public final class SimpleNovelAdjacencyAndChimericAlignmentEvidence {

    private static final NovelAdjacencyAndInferredAltHaptype.Serializer narlSerializer = new NovelAdjacencyAndInferredAltHaptype.Serializer();
    private static final ChimericAlignment.Serializer caSerializer = new ChimericAlignment.Serializer();

    private final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype;
    private final List<ChimericAlignment> alignmentEvidence;

    public NovelAdjacencyAndInferredAltHaptype getNovelAdjacencyReferenceLocations() {
        return novelAdjacencyAndInferredAltHaptype;
    }
    public byte[] getAltHaplotypeSequence() {
        return novelAdjacencyAndInferredAltHaptype.altHaplotypeSequence;
    }
    public List<ChimericAlignment> getAlignmentEvidence() {
        return alignmentEvidence;
    }

    public SimpleNovelAdjacencyAndChimericAlignmentEvidence(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyReferenceLocations,
                                                            final List<ChimericAlignment> alignmentEvidence) {
        this.novelAdjacencyAndInferredAltHaptype = Utils.nonNull( novelAdjacencyReferenceLocations );
        this.alignmentEvidence = Utils.nonNull( alignmentEvidence );
    }

    private SimpleNovelAdjacencyAndChimericAlignmentEvidence(final Kryo kryo, final Input input) {
        novelAdjacencyAndInferredAltHaptype = narlSerializer.read(kryo, input, NovelAdjacencyAndInferredAltHaptype.class);
        final int evidenceCount = input.readInt();
        alignmentEvidence = new ArrayList<>(evidenceCount);
        for (int i = 0; i < evidenceCount; ++i) {
            alignmentEvidence.add(caSerializer.read(kryo, input, ChimericAlignment.class));
        }
    }

    private void serialize(final Kryo kryo, final Output output) {
        narlSerializer.write(kryo, output, novelAdjacencyAndInferredAltHaptype);
        output.writeInt(alignmentEvidence.size());
        alignmentEvidence.forEach(ca -> caSerializer.write(kryo, output, ca));
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleNovelAdjacencyAndChimericAlignmentEvidence> {
        @Override
        public void write(final Kryo kryo, final Output output, final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyReferenceLocations ) {
            novelAdjacencyReferenceLocations.serialize(kryo, output);
        }

        @Override
        public SimpleNovelAdjacencyAndChimericAlignmentEvidence read(final Kryo kryo, final Input input, final Class<SimpleNovelAdjacencyAndChimericAlignmentEvidence> klass ) {
            return new SimpleNovelAdjacencyAndChimericAlignmentEvidence(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SimpleNovelAdjacencyAndChimericAlignmentEvidence that = (SimpleNovelAdjacencyAndChimericAlignmentEvidence) o;

        if (!novelAdjacencyAndInferredAltHaptype.equals(that.novelAdjacencyAndInferredAltHaptype)) return false;
        return alignmentEvidence.equals(that.alignmentEvidence);
    }

    @Override
    public int hashCode() {
        int result = novelAdjacencyAndInferredAltHaptype.hashCode();
        result = 31 * result + alignmentEvidence.hashCode();
        return result;
    }

}
