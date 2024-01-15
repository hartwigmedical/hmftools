package com.hartwig.hmftools.crest;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class CrestApplication
{
    private static final String RNA_SAMPLE = "rna_sample";
    private static final String MIN_TOTAL_READS = "min_total_reads";
    private static final String MIN_RNA_READS = "min_rna_reads";
    private static final String ACCEPTANCE_RATIO = "acceptance_ratio";
    private static final String DO_NOT_WRITE_FILE = "do_not_write_file";

    private static final String DEFAULT_MIN_TOTAL_READS = "10";
    private static final String DEFAULT_MIN_RNA_READS = "1";
    private static final String DEFAULT_ACCEPTANCE_RATIO = "0.9";

    public static void main(String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CrestAlgo crestAlgo = new CrestAlgo(configBuilder.getValue(PURPLE_DIR_CFG),
                parseOutputDir(configBuilder),
                configBuilder.getValue(SAMPLE),
                configBuilder.getValue(RNA_SAMPLE),
                configBuilder.getInteger(MIN_TOTAL_READS),
                configBuilder.getInteger(MIN_RNA_READS),
                configBuilder.getDecimal(ACCEPTANCE_RATIO),
                configBuilder.hasFlag(DO_NOT_WRITE_FILE)
        );

        crestAlgo.run();
    }

    private static void registerConfig(@NotNull final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addConfigItem(RNA_SAMPLE, "ID of RNA sample in vcf to check, e.g. 'COLO829_RNA");
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addFlag(DO_NOT_WRITE_FILE, "Do not write final success or failure file");
        configBuilder.addConfigItem(MIN_TOTAL_READS, false, "Minimum total reads for a SNP to be included", DEFAULT_MIN_TOTAL_READS);
        configBuilder.addConfigItem(MIN_RNA_READS, false, "Minimum allele reads for a SNP to be included", DEFAULT_MIN_RNA_READS);
        configBuilder.addConfigItem(ACCEPTANCE_RATIO, false, "Lower bound on fraction of allele to total reads for check to pass", DEFAULT_ACCEPTANCE_RATIO);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }
}