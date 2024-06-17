package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE_DESC;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.REGIONS_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.checkFileExists;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.from;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.BedFileReader;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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

    public final List<ChrBaseRegion> UnmappableRegions;

    public final Map<String,List<BaseRegion>> TargetRegions;
    public final boolean OnlyTargetRegions;

    // metrics capture config
    public final boolean ExcludeZeroCoverage;
    public final boolean WriteOffTarget;
    public final int HighFragmentOverlapThreshold;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    public final SpecificRegions SpecificChrRegions;

    // debug
    public final List<String> LogReadIds;
    public final boolean PerfDebug;

    private boolean mIsValid;

    private static final String MAP_QUAL_THRESHOLD = "map_qual_threshold";
    private static final String BASE_QUAL_THRESHOLD = "base_qual_threshold";
    private static final String MAX_COVERAGE = "max_coverage";
    private static final String EXCLUDE_ZERO_COVERAGE = "exclude_zero_coverage";
    private static final String ONLY_TARGET = "only_target";

    private static final String OFF_TARGET_FRAG_OVERLAP_THRESHOLD = "off_target_frag_overlap_threshold";
    private static final String WRITE_OFF_TARGET = "write_off_target";

    private static final int DEFAULT_MAP_QUAL_THRESHOLD = 20;
    private static final int DEFAULT_BASE_QUAL_THRESHOLD = 10;
    private static final int DEFAULT_MAX_COVERAGE = 250;

    public MetricsConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        SampleId =  configBuilder.getValue(SAMPLE);
        BamFile =  configBuilder.getValue(BAM_FILE);
        RefGenomeFile =  configBuilder.getValue(REF_GENOME);

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
        }
        else
        {
            OutputDir = pathFromFile(BamFile);
        }

        OutputId =  configBuilder.getValue(OUTPUT_ID);

        if(BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: bam({}) refGenome({}) outputDir({})",
                    BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = configBuilder.hasValue(REF_GENOME_VERSION) ? from(configBuilder) : deriveRefGenomeVersion(BamFile);
        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        BT_LOGGER.info("output({})", OutputDir);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        MapQualityThreshold = configBuilder.getInteger(MAP_QUAL_THRESHOLD);
        BaseQualityThreshold = configBuilder.getInteger(BASE_QUAL_THRESHOLD);
        MaxCoverage = configBuilder.getInteger(MAX_COVERAGE);
        ExcludeZeroCoverage = configBuilder.hasFlag(EXCLUDE_ZERO_COVERAGE);
        WriteOffTarget = configBuilder.hasFlag(WRITE_OFF_TARGET);
        HighFragmentOverlapThreshold = configBuilder.getInteger(OFF_TARGET_FRAG_OVERLAP_THRESHOLD);

        TargetRegions = loadChrBaseRegions(configBuilder.getValue(REGIONS_FILE));
        OnlyTargetRegions = !TargetRegions.isEmpty() && configBuilder.hasFlag(ONLY_TARGET);

        UnmappableRegions = Lists.newArrayList();
        loadUnmappableRegions();

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        mIsValid &= SpecificChrRegions != null;

        LogReadIds = parseLogReadIds(configBuilder);

        Threads = parseThreads(configBuilder);

        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        mIsValid = checkFileExists(BamFile) && checkFileExists(RefGenomeFile);
        return mIsValid;
    }

    private void loadUnmappableRegions()
    {
        String filename = RefGenVersion.is37() ? "/genome_unmappable_regions.37.bed" : "/genome_unmappable_regions.38.bed";

        final InputStream inputStream = MetricsConfig.class.getResourceAsStream(filename);

        try
        {
            UnmappableRegions.addAll(
                    BedFileReader.loadBedFile(new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList())));
        }
        catch(Exception e)
        {
            BT_LOGGER.error("failed to load unmapped regions file({}): {}", filename, e.toString());
            System.exit(1);
        }
    }

    public String formFilename(final String fileType)
    {
        String filename = OutputDir + SampleId + BamMetricsSummary.BAM_METRICS_FILE_ID;

        filename += "." + fileType;

        if(OutputId != null)
            filename += "." + OutputId;

        return filename + TSV_EXTENSION;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, true);
        addRefGenomeVersion(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addPath(BAM_FILE, true, BAM_FILE_DESC);
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);

        configBuilder.addPath(
                REGIONS_FILE, false,
                "TSV or BED file with regions to analyse, expected columns Chromosome,PositionStart,PositionEnd or no headers");

        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addInteger(MAP_QUAL_THRESHOLD, "Map quality threshold", DEFAULT_MAP_QUAL_THRESHOLD);
        configBuilder.addInteger(BASE_QUAL_THRESHOLD, "Base quality threshold", DEFAULT_BASE_QUAL_THRESHOLD);
        configBuilder.addInteger(MAX_COVERAGE, "Max coverage", DEFAULT_MAX_COVERAGE);

        configBuilder.addFlag(ONLY_TARGET, "Only capture metrics within the specific regions file");

        configBuilder.addInteger(
                OFF_TARGET_FRAG_OVERLAP_THRESHOLD,
                "Write regions of high off-target fragment overlap if pile-up above threshold (0=disabled)", 0);

        configBuilder.addFlag(EXCLUDE_ZERO_COVERAGE, "Exclude bases with zero coverage");
        configBuilder.addFlag(WRITE_OFF_TARGET, "Write off-target data");
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    @VisibleForTesting
    public MetricsConfig(int maxCoveage)
    {
        mIsValid = true;

        SampleId = "SAMPLE_ID";
        BamFile = null;
        RefGenomeFile = null;
        RefGenVersion = V37;
        OutputDir = null;
        OutputId = null;

        PartitionSize = DEFAULT_CHR_PARTITION_SIZE;
        MapQualityThreshold = DEFAULT_MAP_QUAL_THRESHOLD;
        BaseQualityThreshold = DEFAULT_BASE_QUAL_THRESHOLD;
        HighFragmentOverlapThreshold = 0;
        MaxCoverage = maxCoveage;
        ExcludeZeroCoverage = false;
        WriteOffTarget = false;

        SpecificChrRegions = new SpecificRegions();
        LogReadIds = Collections.emptyList();
        UnmappableRegions = Collections.emptyList();
        TargetRegions = Maps.newHashMap();
        OnlyTargetRegions = false;

        Threads = 0;
        PerfDebug = false;
    }
}
