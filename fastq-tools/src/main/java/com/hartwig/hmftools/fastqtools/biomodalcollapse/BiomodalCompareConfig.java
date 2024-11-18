package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BiomodalCompareConfig
{
    public final String RefBam;
    public final String NewBam;
    public final String OutputTsv;

    private static final String REF_BAM = "ref_bam";
    private static final String NEW_BAM = "new_bam";
    private static final String OUTPUT_TSV = "output_tsv";

    public BiomodalCompareConfig(final ConfigBuilder configBuilder)
    {
        RefBam = configBuilder.getValue(REF_BAM);
        NewBam = configBuilder.getValue(NEW_BAM);
        OutputTsv = configBuilder.getValue(OUTPUT_TSV);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REF_BAM, true, "Path to the reference bam file");
        configBuilder.addConfigItem(NEW_BAM, true, "Path to the new bam file");
        configBuilder.addConfigItem(OUTPUT_TSV, true, "Pather to output TSV file");

        addLoggingOptions(configBuilder);
    }
}
