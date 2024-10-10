package com.hartwig.hmftools.bamtools.fastqbiomodalcollapse;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FastqBiomodalCollapseStatsConfig
{
    public final String SampleId;
    public final String Fastq1Path;
    public final String Fastq2Path;
    public final String OutputDir;

    private static final String SAMPLE_ID = "sample";
    private static final String FASTQ1_PATH = "fastq1";
    private static final String FASTQ2_PATH = "fastq2";
    private static final String OUTPUT_DIR = "output_dir";

    public FastqBiomodalCollapseStatsConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE_ID);
        Fastq1Path = configBuilder.getValue(FASTQ1_PATH);
        Fastq2Path = configBuilder.getValue(FASTQ2_PATH);
        OutputDir = configBuilder.getValue(OUTPUT_DIR);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE_ID, true, "ID of the sample");
        configBuilder.addConfigItem(FASTQ1_PATH, true, "Path to the fastq file of first in pair reads");
        configBuilder.addConfigItem(FASTQ2_PATH, true, "Path to the fastq file of second in pair reads");
        configBuilder.addConfigItem(OUTPUT_DIR, true, "Directory for output files");
    }
}
