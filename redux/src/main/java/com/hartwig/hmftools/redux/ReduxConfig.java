package com.hartwig.hmftools.redux;

import static java.lang.Math.ceil;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig.JITTER_MSI_SITES_FILE;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig.JITTER_MSI_SITES_FILE_DESC;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.UnmappedRegions.UNMAP_REGIONS_FILE;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SEQUENCING_TYPE_CFG;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_LOG_TIME;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_LOG_TIME_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.redux.common.Constants.FILE_ID;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MIN_HIGH_DEPTH;
import static com.hartwig.hmftools.redux.write.ReadOutput.NONE;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.region.UnmappedRegions;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.redux.common.FilterReadsType;
import com.hartwig.hmftools.redux.unmap.ReadUnmapper;
import com.hartwig.hmftools.redux.umi.UmiConfig;
import com.hartwig.hmftools.redux.write.ReadOutput;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.ValidationStringency;

public class ReduxConfig
{
    public final String SampleId;
    public final List<String> BamFiles;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeInterface RefGenome;
    public final SequencingType Sequencing;

    public final ValidationStringency BamStringency;

    // UMI group config
    public final UmiConfig UMIs;
    public final boolean FormConsensus;
    public final boolean JitterMsiOnly;

    public final ReadUnmapper UnmapRegions;

    public final String OutputBam;
    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteBam;
    public final boolean MultiBam;
    public final boolean WriteStats;

    public final int Threads;
    public final int PartitionThreadRatio;

    public final String BamToolPath;

    // debug
    public final boolean KeepInterimBams;
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final FilterReadsType SpecificRegionsFilterType;
    public final ReadOutput LogReadType;
    public final double PerfDebugTime;
    public final boolean RunChecks;
    public final boolean DropDuplicates;
    public final boolean LogFinalCache;
    public final int WriteReadBaseLength;

    public final String JitterMsiFile;
    public final double JitterMsiMaxSitePercContribution;

    private boolean mIsValid;
    private int mReadLength;

    public static final Logger RD_LOGGER = LogManager.getLogger(ReduxConfig.class);
    public static final String APP_NAME = "Redux";

    // config strings
    private static final String INPUT_BAM = "input_bam";
    private static final String BAM_FILE = "bam_file";
    private static final String OUTPUT_BAM = "output_bam";
    private static final String READ_OUTPUTS = "read_output";
    private static final String FORM_CONSENSUS = "form_consensus";
    private static final String READ_LENGTH = "read_length";

    private static final String WRITE_STATS = "write_stats";
    private static final String DROP_DUPLICATES = "drop_duplicates";
    private static final String JITTER_MSI_ONLY = "jitter_msi_only";
    private static final String PARTIION_THREAD_RATIO = "partition_ratio";

    private static final String JITTER_MSI_MAX_SINGLE_SITE_ALT_CONTRIBUTION = "jitter_max_site_alt_contribution";

    // dev and options
    public static final String KEEP_INTERIM_BAMS = "keep_interim_bams";
    private static final String NO_WRITE_BAM = "no_write_bam";
    private static final String RUN_CHECKS = "run_checks";
    private static final String LOG_FINAL_CACHE = "log_final_cache";
    private static final String SPECIFIC_REGION_FILTER_TYPE = "specific_region_filter";
    private static final String WRITE_READ_BASE_LENGTH = "write_read_base_length";

    public ReduxConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;
        SampleId = configBuilder.getValue(SAMPLE);

        String bamFiles = configBuilder.hasValue(INPUT_BAM) ? configBuilder.getValue(INPUT_BAM) : configBuilder.getValue(BAM_FILE);
        BamFiles = Arrays.stream(bamFiles.split(CONFIG_FILE_DELIM, -1)).collect(Collectors.toList());

        if(BamFiles.isEmpty())
        {
            RD_LOGGER.error("no BAM files configured");
            System.exit(1);
        }

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);

        Sequencing = SequencingType.valueOf(configBuilder.getValue(SEQUENCING_TYPE_CFG));

        OutputBam = configBuilder.getValue(OUTPUT_BAM);

        if(OutputBam != null && BamFiles.stream().anyMatch(x -> x.equals(OutputBam)))
        {
            RD_LOGGER.error("output BAM({}) matches input BAM filename", OutputBam);
            System.exit(1);
        }

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
        }
        else if(OutputBam != null)
        {
            OutputDir = pathFromFile(OutputBam);
        }
        else
        {
            OutputDir = pathFromFile(BamFiles.get(0));
        }

        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(SampleId == null || OutputDir == null || RefGenomeFile == null)
        {
            RD_LOGGER.error("missing config: sample({}) refGenome({}) outputDir({})",
                    SampleId != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        // MD_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        RD_LOGGER.info("output({})", OutputDir);

        BamStringency = BamUtils.validationStringency(configBuilder);

        mReadLength = configBuilder.getInteger(READ_LENGTH);

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        UMIs = UmiConfig.from(configBuilder);

        FormConsensus = UMIs.Enabled || configBuilder.hasFlag(FORM_CONSENSUS);

        if(configBuilder.hasValue(UNMAP_REGIONS_FILE))
        {
            UnmapRegions = new ReadUnmapper(configBuilder.getValue(UNMAP_REGIONS_FILE));
        }
        else
        {
            Map<String, List<HighDepthRegion>> unmappedMap = Maps.newHashMap();

            ChrBaseRegion excludedRegion = ExcludedRegions.getPolyGRegion(RefGenVersion);
            unmappedMap.put(excludedRegion.Chromosome, Lists.newArrayList(HighDepthRegion.from(excludedRegion, UNMAP_MIN_HIGH_DEPTH)));

            UnmapRegions = new ReadUnmapper(unmappedMap);
        }

        String duplicateLogic = UMIs.Enabled ? "UMIs" : (FormConsensus ? "consensus" : "max base-qual");
        RD_LOGGER.info("duplicate logic: {}", duplicateLogic);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;

        SpecificRegionsFilterType = SpecificChrRegions.hasFilters() ?
                FilterReadsType.valueOf(configBuilder.getValue(SPECIFIC_REGION_FILTER_TYPE, FilterReadsType.READ.toString())) :
                FilterReadsType.NONE;

        Threads = parseThreads(configBuilder);
        PartitionThreadRatio = configBuilder.getInteger(PARTIION_THREAD_RATIO);

        LogReadType = ReadOutput.valueOf(configBuilder.getValue(READ_OUTPUTS, NONE.toString()));

        JitterMsiOnly = configBuilder.hasFlag(JITTER_MSI_ONLY);
        JitterMsiFile = configBuilder.getValue(JITTER_MSI_SITES_FILE);
        JitterMsiMaxSitePercContribution = configBuilder.getDecimal(JITTER_MSI_MAX_SINGLE_SITE_ALT_CONTRIBUTION);

        WriteBam = !configBuilder.hasFlag(NO_WRITE_BAM) && !JitterMsiOnly;
        MultiBam = WriteBam && Threads > 1; // now on automatically
        KeepInterimBams = configBuilder.hasFlag(KEEP_INTERIM_BAMS);

        LogReadIds = parseLogReadIds(configBuilder);

        WriteStats = configBuilder.hasFlag(WRITE_STATS);
        PerfDebugTime = configBuilder.getDecimal(PERF_LOG_TIME);
        RunChecks = configBuilder.hasFlag(RUN_CHECKS);
        LogFinalCache = configBuilder.hasFlag(LOG_FINAL_CACHE);
        DropDuplicates = configBuilder.hasFlag(DROP_DUPLICATES);
        WriteReadBaseLength = configBuilder.getInteger(WRITE_READ_BASE_LENGTH);

        if(RunChecks)
        {
            RD_LOGGER.info("running debug options: read-checks({})", RunChecks);
        }
    }

    public boolean isValid() { return mIsValid; }

    public int readLength() { return mReadLength; }
    public void setReadLength(int readLength)
    {
        mReadLength = readLength;
    }

    public boolean perfDebug() { return PerfDebugTime > 0; }

    public String formFilename(final String fileType)
    {
        String filename = OutputDir + SampleId + "." + FILE_ID;

        filename += "." + fileType;

        if(OutputId != null)
            filename += "." + OutputId;

        filename += TSV_EXTENSION;

        return filename;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, false, "BAM filename, deprecated, use 'input_bam' instead");
        configBuilder.addPaths(INPUT_BAM, false, "BAM file path, separated by ',' if multiple");
        configBuilder.addConfigItem(OUTPUT_BAM, false, "Output BAM filename");
        addRefGenomeConfig(configBuilder, true);
        SequencingType.registerConfig(configBuilder);
        configBuilder.addInteger(READ_LENGTH, "Read length, otherwise will sample from BAM", 0);

        configBuilder.addInteger(PARTIION_THREAD_RATIO, "Partitions per thread, impacts BAM-writing performance", 2);

        configBuilder.addConfigItem(
                READ_OUTPUTS, false, format("Write reads: %s", ReadOutput.valuesStr()), NONE.toString());

        UnmappedRegions.registerConfig(configBuilder);

        configBuilder.addFlag(NO_WRITE_BAM, "BAM not written, producing only TSV reads and/or statistics");
        configBuilder.addFlag(KEEP_INTERIM_BAMS, "Do no delete per-thread BAMs");

        BamToolName.addConfig(configBuilder);

        configBuilder.addFlag(FORM_CONSENSUS, "Form consensus reads from duplicate groups without UMIs");
        configBuilder.addFlag(WRITE_STATS, "Write duplicate and UMI-group stats");
        configBuilder.addFlag(DROP_DUPLICATES, "Drop duplicates from BAM");
        configBuilder.addFlag(JITTER_MSI_ONLY, "Jitter MSi output only, no duplicate processing");
        addValidationStringencyOption(configBuilder);
        UmiConfig.addConfig(configBuilder);

        configBuilder.addPath(JITTER_MSI_SITES_FILE, false, JITTER_MSI_SITES_FILE_DESC);

        configBuilder.addDecimal(
                JITTER_MSI_MAX_SINGLE_SITE_ALT_CONTRIBUTION, "Jitter MIS max single alt site perc contribute",
                DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION);

        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, PERF_LOG_TIME_DESC, 0);
        configBuilder.addFlag(RUN_CHECKS, "Run duplicate mismatch checks");
        configBuilder.addFlag(LOG_FINAL_CACHE, "Log cached fragments on completion");
        configBuilder.addConfigItem(SPECIFIC_REGION_FILTER_TYPE, "Used with specific regions, to filter mates or supps");

        configBuilder.addInteger(WRITE_READ_BASE_LENGTH, "Number of read bases to write with read data", 0);
    }

    @VisibleForTesting
    public ReduxConfig(
            final RefGenomeInterface refGenome, boolean umiEnabled, boolean duplexUmi, boolean formConsensus, final ReadUnmapper readUnmapper)
    {
        mIsValid = true;
        SampleId = "";
        BamFiles = Lists.newArrayList();
        RefGenomeFile = null;
        OutputBam = null;
        OutputDir = null;
        OutputId = "";
        RefGenVersion = V37;
        RefGenome = refGenome;
        Sequencing = ILLUMINA;

        BamStringency = ValidationStringency.STRICT;
        mReadLength = DEFAULT_READ_LENGTH;

        UMIs = new UmiConfig(umiEnabled, duplexUmi, String.valueOf(DEFAULT_DUPLEX_UMI_DELIM), false);
        FormConsensus = umiEnabled || formConsensus;

        SpecificChrRegions = new SpecificRegions();
        SpecificRegionsFilterType = FilterReadsType.MATE_AND_SUPP;

        BamToolPath = null;

        UnmapRegions = readUnmapper;

        JitterMsiFile = null;
        JitterMsiMaxSitePercContribution = DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION;
        JitterMsiOnly = false;

        WriteBam = false;
        MultiBam = false;
        KeepInterimBams = false;
        LogReadType = NONE;

        LogReadIds = Lists.newArrayList();
        Threads = 0;
        PartitionThreadRatio = 1;
        PerfDebugTime = 0;
        RunChecks = true;
        WriteStats = false;
        LogFinalCache = true;
        DropDuplicates = false;
        WriteReadBaseLength = 0;
    }
}
