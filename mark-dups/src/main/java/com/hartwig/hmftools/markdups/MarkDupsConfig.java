package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_PARTITION_SIZE;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_POS_BUFFER_SIZE;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.markdups.common.Constants.UNMAP_MIN_HIGH_DEPTH;
import static com.hartwig.hmftools.markdups.write.ReadOutput.NONE;

import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.samtools.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.markdups.common.FilterReadsType;
import com.hartwig.hmftools.markdups.common.HighDepthRegion;
import com.hartwig.hmftools.markdups.common.ReadUnmapper;
import com.hartwig.hmftools.markdups.consensus.GroupIdGenerator;
import com.hartwig.hmftools.markdups.umi.UmiConfig;
import com.hartwig.hmftools.markdups.write.ReadOutput;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.ValidationStringency;

public class MarkDupsConfig
{
    public final String SampleId;
    public final List<String> BamFiles;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeInterface RefGenome;

    public final int PartitionSize;
    public final int BufferSize;
    public final ValidationStringency BamStringency;

    // UMI group config
    public final UmiConfig UMIs;
    public final boolean FormConsensus;
    public final GroupIdGenerator IdGenerator;

    public final ReadUnmapper UnmapRegions;

    public final String OutputBam;
    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteBam;
    public final boolean MultiBam;
    public final boolean WriteStats;

    public final boolean NoMateCigar;
    public final int Threads;

    public final String SamToolsPath;
    public final String SambambaPath;

    // debug
    public final boolean KeepInterimBams;
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final FilterReadsType SpecificRegionsFilterType;
    public final ReadOutput LogReadType;
    public final boolean PerfDebug;
    public final boolean RunChecks;
    public final boolean LogFinalCache;

    private boolean mIsValid;
    private int mReadLength;

    public static final Logger MD_LOGGER = LogManager.getLogger(MarkDupsConfig.class);
    public static final String APP_NAME = "MarkDups";

    // config strings
    private  static final String INPUT_BAM = "input_bam";
    private  static final String BAM_FILE = "bam_file";
    private  static final String OUTPUT_BAM = "output_bam";
    private static final String PARTITION_SIZE = "partition_size";
    private static final String BUFFER_SIZE = "buffer_size";
    private static final String READ_OUTPUTS = "read_output";
    private static final String NO_MATE_CIGAR = "no_mate_cigar";
    private static final String FORM_CONSENSUS = "form_consensus";
    private static final String READ_LENGTH = "read_length";

    private static final String MULTI_BAM = "multi_bam";
    private static final String SAMTOOLS_PATH = "samtools";
    private static final String SAMBAMBA_PATH = "sambamba";
    private static final String UNMAP_REGIONS = "unmap_regions";
    private static final String WRITE_STATS = "write_stats";

    // debug
    private static final String KEEP_INTERIM_BAMS = "keep_interim_bams";
    private static final String NO_WRITE_BAM = "no_write_bam";
    private static final String RUN_CHECKS = "run_checks";
    private static final String LOG_FINAL_CACHE = "log_final_cache";
    private static final String SPECIFIC_REGION_FILTER_TYPE = "specific_region_filter";

    public MarkDupsConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;
        SampleId = configBuilder.getValue(SAMPLE);

        String bamFiles = configBuilder.hasValue(INPUT_BAM) ? configBuilder.getValue(INPUT_BAM) : configBuilder.getValue(BAM_FILE);
        BamFiles = Arrays.stream(bamFiles.split(CONFIG_FILE_DELIM, -1)).collect(Collectors.toList());

        if(BamFiles.isEmpty())
        {
            MD_LOGGER.error("no BAM files configured");
            System.exit(1);
        }

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);

        OutputBam = configBuilder.getValue(OUTPUT_BAM);

        if(OutputBam != null && BamFiles.stream().anyMatch(x -> x.equals(OutputBam)))
        {
            MD_LOGGER.error("output BAM({}) matches input BAM filename", OutputBam);
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

        if(SampleId == null|| OutputDir == null || RefGenomeFile == null)
        {
            MD_LOGGER.error("missing config: sample({}) refGenome({}) outputDir({})",
                    SampleId != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        // MD_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        MD_LOGGER.info("output({})", OutputDir);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        BufferSize = configBuilder.getInteger(BUFFER_SIZE);
        BamStringency = BamUtils.validationStringency(configBuilder);

        mReadLength = configBuilder.getInteger(READ_LENGTH);

        SambambaPath = configBuilder.getValue(SAMBAMBA_PATH);
        SamToolsPath = configBuilder.getValue(SAMTOOLS_PATH);

        NoMateCigar = configBuilder.hasFlag(NO_MATE_CIGAR);
        UMIs = UmiConfig.from(configBuilder);
        FormConsensus = !UMIs.Enabled && !NoMateCigar && configBuilder.hasFlag(FORM_CONSENSUS);
        IdGenerator = new GroupIdGenerator();

        if(configBuilder.hasValue(UNMAP_REGIONS))
        {
            UnmapRegions = new ReadUnmapper(configBuilder.getValue(UNMAP_REGIONS));
        }
        else
        {
            Map<String, List<HighDepthRegion>> unmappedMap = Maps.newHashMap();

            ChrBaseRegion excludedRegion = ExcludedRegions.getPolyGRegion(RefGenVersion);
            unmappedMap.put(excludedRegion.Chromosome, Lists.newArrayList(HighDepthRegion.from(excludedRegion, UNMAP_MIN_HIGH_DEPTH)));

            UnmapRegions = new ReadUnmapper(unmappedMap);
        }

        String duplicateLogic = UMIs.Enabled ? "UMIs" : (FormConsensus ? "consensus" : "max base-qual");
        MD_LOGGER.info("duplicate logic: {}", duplicateLogic);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;

        SpecificRegionsFilterType = SpecificChrRegions.hasFilters() ?
                FilterReadsType.valueOf(configBuilder.getValue(SPECIFIC_REGION_FILTER_TYPE, FilterReadsType.READ.toString())) :
                FilterReadsType.NONE;

        Threads = parseThreads(configBuilder);

        LogReadType = ReadOutput.valueOf(configBuilder.getValue(READ_OUTPUTS, NONE.toString()));

        WriteBam = !configBuilder.hasFlag(NO_WRITE_BAM);
        MultiBam = WriteBam && Threads > 1 && configBuilder.hasFlag(MULTI_BAM);
        KeepInterimBams = configBuilder.hasFlag(KEEP_INTERIM_BAMS);

        LogReadIds = parseLogReadIds(configBuilder);

        WriteStats = configBuilder.hasFlag(WRITE_STATS);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
        RunChecks = configBuilder.hasFlag(RUN_CHECKS);
        LogFinalCache = configBuilder.hasFlag(LOG_FINAL_CACHE);

        if(RunChecks)
        {
            MD_LOGGER.info("running debug options: read-checks({})", RunChecks);
        }
    }

    public boolean isValid() { return mIsValid; }

    public int readLength() { return mReadLength; }

    public void setReadLength(int readLength)
    {
        mReadLength = readLength;
    }

    public String formFilename(final String fileType)
    {
        String filename = OutputDir + SampleId;

        filename += "." + fileType;

        if(OutputId != null)
            filename += "." + OutputId;

        filename += TSV_EXTENSION;

        return filename;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, false, "BAM filename, deprecated, use 'input_bam' instead");
        configBuilder.addPaths(INPUT_BAM, false, "BAM file path, separated by ',' if multiple");
        configBuilder.addConfigItem(OUTPUT_BAM, false, "Output BAM filename");
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_PARTITION_SIZE);
        configBuilder.addInteger(BUFFER_SIZE, "Read buffer size", DEFAULT_POS_BUFFER_SIZE);
        configBuilder.addInteger(READ_LENGTH, "Read length, otherwise will sample from BAM", 0);

        configBuilder.addConfigItem(
                READ_OUTPUTS, false, format("Write reads: %s", ReadOutput.valuesStr()), NONE.toString());

        configBuilder.addPath(UNMAP_REGIONS, false, "Unmap reads within these regions");

        configBuilder.addFlag(NO_WRITE_BAM, "BAM not written, producing only TSV reads and/or statistics");
        configBuilder.addFlag(MULTI_BAM, "Write temporary BAMs with multi-threading");
        configBuilder.addFlag(KEEP_INTERIM_BAMS, "Do no delete per-thread BAMs");
        configBuilder.addPath(SAMTOOLS_PATH, false, "Path to samtools for sort");
        configBuilder.addPath(SAMBAMBA_PATH, false, "Path to sambamba for merge");

        configBuilder.addFlag(FORM_CONSENSUS, "Form consensus reads from duplicate groups without UMIs");
        configBuilder.addFlag(NO_MATE_CIGAR, "Mate CIGAR not set by aligner, make no attempt to use it");
        configBuilder.addFlag(WRITE_STATS, "Write duplicate and UMI-group stats");
        addValidationStringencyOption(configBuilder);
        UmiConfig.addConfig(configBuilder);

        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addFlag(RUN_CHECKS, "Run duplicate mismatch checks");
        configBuilder.addFlag(LOG_FINAL_CACHE, "Log cached fragments on completion");
        configBuilder.addConfigItem(SPECIFIC_REGION_FILTER_TYPE, "Used with specific regions, to filter mates or supps");
    }

    public MarkDupsConfig(
            int partitionSize, int bufferSize, final RefGenomeInterface refGenome, boolean umiEnabled, boolean duplexUmi, boolean formConsensus)
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

        PartitionSize = partitionSize;
        BufferSize = bufferSize;
        BamStringency = ValidationStringency.STRICT;
        mReadLength = DEFAULT_READ_LENGTH;

        UMIs = new UmiConfig(umiEnabled, duplexUmi, String.valueOf(DEFAULT_DUPLEX_UMI_DELIM), false);
        FormConsensus = formConsensus;
        NoMateCigar = false;
        IdGenerator = new GroupIdGenerator();

        SpecificChrRegions = new SpecificRegions();
        SpecificRegionsFilterType = FilterReadsType.MATE_AND_SUPP;

        SamToolsPath = null;
        SambambaPath = null;

        UnmapRegions = new ReadUnmapper(Maps.newHashMap());

        WriteBam = false;
        MultiBam = false;
        KeepInterimBams = false;
        LogReadType = NONE;

        LogReadIds = Lists.newArrayList();
        Threads = 0;
        PerfDebug = false;
        RunChecks = true;
        WriteStats = false;
        LogFinalCache = true;
    }
}
