package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import java.io.Closeable;
import java.util.List;

public interface AnnotatedGenomicRegionWriter extends Closeable {

    /**
     * Write only the header (and any SAMFileHeader or comments)
     *
     * @param SAMFileHeader
     * @param comments
     * @param annotations
     * @param contigColumnName
     * @param startColumnName
     * @param endColumnName
     */
    void writeHeader(final String SAMFileHeader, final List<String> comments, final List<String> annotations,
                     final String contigColumnName, final String startColumnName, final String endColumnName);

    /**
     * attempt to close the VCF file
     */
    @Override
    void close() ;

    /** Write one region to the file.
     *
     * @param simpleAnnotatedGenomicRegion
     */
    void add(final SimpleAnnotatedGenomicRegion simpleAnnotatedGenomicRegion);
}
