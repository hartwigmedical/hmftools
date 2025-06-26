package com.hartwig.hmftools.chord;

import com.google.common.io.Resources;

public class ChordTestDataPaths
{
    public static final String INPUT_VCF_DIR = Resources.getResource("vcf/").getPath();
    public static final String TMP_OUTPUT_DIR = System.getProperty("java.io.tmpdir") + "/chord_output/";

    public static final String EMPTY_SAMPLE = "EMPTY_SAMPLE";
    public static final String MINIMAL_SAMPLE = "MINIMAL_SAMPLE";
    public static final String NON_STANDARD_NUC_GENOME_SAMPLE = "NON_STANDARD_NUC_GENOME_SAMPLE";

    public static final String DUMMY_GENOME_FASTA = Resources.getResource("fasta/dummy_genome.fasta").getPath();
    public static final String NON_STANDARD_NUC_GENOME_FASTA = Resources.getResource("fasta/non_standard_nuc_genome.fasta").getPath();

}
