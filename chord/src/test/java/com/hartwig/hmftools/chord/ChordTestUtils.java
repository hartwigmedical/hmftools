package com.hartwig.hmftools.chord;

import com.google.common.io.Resources;

public class ChordTestUtils
{
    public static final String EMPTY_SAMPLE = "EMPTY_SAMPLE";
    public static final String MINIMAL_SAMPLE = "MINIMAL_SAMPLE";

    public static final String INPUT_VCF_DIR = Resources.getResource("vcf/").getPath();
    public static final String TMP_OUTPUT_DIR = System.getProperty("java.io.tmpdir") + "/chord_output/";

    public static final String HUMAN_GENOME_FASTA = "/Users/lnguyen/Hartwig/hartwigmedical/resources/ref_genomes/GRCh37/Homo_sapiens.GRCh37.GATK.illumina.fasta";
}
