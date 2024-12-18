package com.hartwig.hmftools.bamtools.copyfastqtags;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class CopyFastqTagsConfig
{
    public final String OutputBamFile;
    public final String BamFile;
    public final String FastqFile;

    private static final String OUTPUT_BAM_FILE = "output";
    private static final String BAM_FILE = "bam_file";
    private static final String FASTQ_FILE = "fastq_file";

    public CopyFastqTagsConfig(final ConfigBuilder configBuilder)
    {
        OutputBamFile = configBuilder.getValue(OUTPUT_BAM_FILE);
        BamFile = configBuilder.getValue(BAM_FILE);
        FastqFile = configBuilder.getValue(FASTQ_FILE);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_BAM_FILE, true, "Output bam file");
        configBuilder.addPath(BAM_FILE, true, "Input bam file");
        configBuilder.addPath(FASTQ_FILE, true, "Input fastq file");

        addLoggingOptions(configBuilder);
    }
}
