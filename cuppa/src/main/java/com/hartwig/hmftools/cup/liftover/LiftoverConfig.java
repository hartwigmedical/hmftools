package com.hartwig.hmftools.cup.liftover;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class LiftoverConfig
{
    public final String OutputDir;
    public final String SampleVcfDir;
    public final int Threads;
    public final boolean ApplyFilters;
    public final boolean KeepExisting;

    public static final String SAMPLE_VCF_DIR = "sample_vcf_dir";
    public static final String APPLY_FILTERS = "apply_filters";
    public static final String KEEP_EXISTING = "keep_existing";

    public LiftoverConfig(final ConfigBuilder configBuilder)
    {
        OutputDir = parseOutputDir(configBuilder);
        SampleVcfDir = configBuilder.getValue(SAMPLE_VCF_DIR);
        Threads = parseThreads(configBuilder);
        ApplyFilters = configBuilder.hasFlag(APPLY_FILTERS);
        KeepExisting = configBuilder.hasFlag(KEEP_EXISTING);
    }

    public static void addOptions(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, false);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
        GenomeLiftoverCache.addConfig(configBuilder);

        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(SAMPLE_VCF_DIR, true, "Path to sample VCF(s)");
        configBuilder.addFlag(APPLY_FILTERS, "Ignore MNVs and non-PASS variants");
        configBuilder.addFlag(KEEP_EXISTING, "Do not overwrite existing output sample files");
    }
}
