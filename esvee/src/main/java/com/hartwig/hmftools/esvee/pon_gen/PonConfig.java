package com.hartwig.hmftools.esvee.pon_gen;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class PonConfig
{
    public final String OutputFilenameSuffix;
    public final String VcfPath;
    public final String OutputDir;

    public final int MinSamples;
    public final int PositionBuffer;
    public final boolean WriteBedFiles;

    public final int Threads;

    public final SpecificRegions SpecificChrRegions;

    private static final String VCF_PATH = "vcf_path";
    private static final String OUTPUT_PON_SUFFIX = "output_pon_suffix";
    private static final String MIN_SAMPLES = "min_samples";
    private static final String POSITION_BUFFER = "position_buffer";
    private static final String WRITE_BED_FILES = "write_bed";

    // constants
    private static final int DEFAULT_MIN_SAMPLES = 3;

    public PonConfig(final ConfigBuilder configBuilder)
    {
        OutputDir = parseOutputDir(configBuilder);
        VcfPath = configBuilder.getValue(VCF_PATH);
        OutputFilenameSuffix = configBuilder.getValue(OUTPUT_PON_SUFFIX);

        MinSamples = configBuilder.getInteger(MIN_SAMPLES);
        PositionBuffer = configBuilder.getInteger(POSITION_BUFFER);

        WriteBedFiles = configBuilder.hasFlag(WRITE_BED_FILES);

        Threads = parseThreads(configBuilder);

        SpecificChrRegions = SpecificRegions.from(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, true);
        configBuilder.addConfigItem(VCF_PATH, true, "VCF path for samples");
        configBuilder.addConfigItem(OUTPUT_PON_SUFFIX, false, "Output PON filename");

        configBuilder.addInteger(MIN_SAMPLES, "Min samples for variant to be included in PON", DEFAULT_MIN_SAMPLES);
        configBuilder.addInteger(POSITION_BUFFER, "Merge PON entries if positions are within buffer", 0);

        configBuilder.addFlag(WRITE_BED_FILES, "Write output in BED format instead of TSV");

        addThreadOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addOutputDir(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
