package com.hartwig.hmftools.common.variant.vcfloader;

final class VCFTestConstants {

    private VCFTestConstants() {
    }

    static final String COMMENT_LINE = "## -> This is a comment line in VCF!";
    static final String HEADER_LINE = "#CHROM -> This line is a header line in VCF!";
    static final String PASS_DATA_LINE_1 = "Col1 \t Col2 \t Col3 \t Col4 \t Col5 \t Col6 \tPASS\t-> Passing!";
    static final String PASS_DATA_LINE_2 = "Col1 \t Col2 \t Col3 \t Col4 \t Col5 \t Col6 \t.\t-> Passing!";
    static final String FILTERED_DATA_LINE = "Col1 \t Col2 \t Col3 \t Col4 \t Col5 \t Col6 \t???\t-> Filtered!";
}
