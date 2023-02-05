package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.checkFileExists;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.CommonUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MetricsConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int MapQualityThreshold;
    public final int BaseQualityThreshold;
    public final int MaxCoverage;

    public final int PartitionSize;

    // metrics capture config
    public final boolean ExcludeZeroCoverage;
    public final boolean WriteOldStyle;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean PerfDebug;

    private boolean mIsValid;

    public static final Logger BT_LOGGER = LogManager.getLogger(MetricsConfig.class);

    // config strings
    public static final String SAMPLE = "sample";
    public static final String BAM_FILE = "bam_file";
    public static final String PARTITION_SIZE = "partition_size";

    private static final String MAP_QUAL_THRESHOLD = "map_qual_threshold";
    private static final String BASE_QUAL_THRESHOLD = "base_qual_threshold";
    private static final String MAX_COVERAGE = "max_coverage";
    private static final String EXCLUDE_ZERO_COVERAGE = "exclude_zero_coverage";
    private static final String WRITE_OLD_STYLE = "write_old_style";

    public static final String LOG_READ_IDS = "log_read_ids";
    public static final String PERF_DEBUG = "perf_debug";

    public static final String ITEM_DELIM = ";";

    private static final int DEFAULT_MAP_QUAL_THRESHOLD = 20;
    private static final int DEFAULT_BASE_QUAL_THRESHOLD = 10;
    private static final int DEFAULT_MAX_COVERAGE = 250;

    public MetricsConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        if(SampleId == null || BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    SampleId != null, BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        BT_LOGGER.info("output({})", OutputDir);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));
        MapQualityThreshold = Integer.parseInt(cmd.getOptionValue(MAP_QUAL_THRESHOLD, String.valueOf(DEFAULT_MAP_QUAL_THRESHOLD)));
        BaseQualityThreshold = Integer.parseInt(cmd.getOptionValue(BASE_QUAL_THRESHOLD, String.valueOf(DEFAULT_BASE_QUAL_THRESHOLD)));
        MaxCoverage = Integer.parseInt(cmd.getOptionValue(MAX_COVERAGE, String.valueOf(DEFAULT_MAX_COVERAGE)));
        ExcludeZeroCoverage = cmd.hasOption(EXCLUDE_ZERO_COVERAGE);
        WriteOldStyle = cmd.hasOption(WRITE_OLD_STYLE);

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();

        Threads = parseThreads(cmd);

        PerfDebug = cmd.hasOption(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        mIsValid = checkFileExists(BamFile) && checkFileExists(RefGenomeFile);
        return mIsValid;
    }

    public String formFilename(final String fileType) { return CommonUtils.formFilename(SampleId, OutputDir, OutputId, fileType); }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "BAM file location");
        addRefGenomeConfig(options);;
        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_CHR_PARTITION_SIZE);
        options.addOption(MAP_QUAL_THRESHOLD, true, "Map quality threshold, default: " + DEFAULT_MAP_QUAL_THRESHOLD);
        options.addOption(BASE_QUAL_THRESHOLD, true, "Base quality threshold, default: " + DEFAULT_BASE_QUAL_THRESHOLD);
        options.addOption(MAX_COVERAGE, true, "Max coverage, default: " + DEFAULT_MAX_COVERAGE);
        options.addOption(EXCLUDE_ZERO_COVERAGE, false, "Exclude bases with zero coverage");
        options.addOption(WRITE_OLD_STYLE, false, "Write data in same format as Picard CollectWgsMetrics");
        addThreadOptions(options);

        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(PERF_DEBUG, false, "Detailed performance tracking and logging");

        return options;
    }
}
