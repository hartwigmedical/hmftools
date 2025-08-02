package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.mappability.UnmappedRegions.UNMAP_REGIONS_FILE;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SEQUENCING_TYPE_CFG;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;
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
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.redux.common.Constants.FILE_ID;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MIN_HIGH_DEPTH;
import static com.hartwig.hmftools.redux.write.ReadOutput.NONE;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.redux.bqr.BqrConfig;
import com.hartwig.hmftools.redux.jitter.JitterAnalyserConfig;
import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.mappability.UnmappedRegions;
import com.hartwig.hmftools.common.mappability.UnmappingRegion;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.redux.common.DuplicateGroupCollapseConfig;
import com.hartwig.hmftools.redux.common.FilterReadsType;
import com.hartwig.hmftools.redux.umi.UmiConfig;
import com.hartwig.hmftools.redux.unmap.ReadChecker;
import com.hartwig.hmftools.redux.unmap.ReadUnmapper;
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
    public final JitterAnalyserConfig JitterConfig;

    // global for convenience
    public static SequencingType SEQUENCING_TYPE = ILLUMINA;

    public final ValidationStringency BamStringency;

    public final DuplicateGroupCollapseConfig DuplicateGroupCollapse;

    // UMI group config
    public final UmiConfig UMIs;
    public final boolean FormConsensus;
    public final boolean JitterMsiOnly;

    public final BqrConfig BQR;

    public final ReadUnmapper UnmapRegions;
    public final boolean SkipUnmapping; // to skip unmapping in-built excluded regions

    public final String OutputBam;
    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteBam;
    public final boolean MultiBam;
    public final boolean WriteStats;

    public final int Threads;
    public final int PartitionThreadRatio;

    public final String BamToolPath;
    public final boolean ParallelConcatenation;

    // debug
    public final boolean KeepInterimBams;
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final FilterReadsType SpecificRegionsFilterType;
    public final ReadOutput LogReadType;
    public final double PerfDebugTime;
    public final boolean FailOnMissingSuppMateCigar;
    public final boolean RunChecks;
    public final boolean DropDuplicates;
    public final boolean SkipFullyUnmappedReads;
    public final int WriteReadBaseLength;

    private boolean mIsValid;
    private int mReadLength;
    private final ReadChecker mReadChecker;

    public static final Logger RD_LOGGER = LogManager.getLogger(ReduxConfig.class);
    public static final String APP_NAME = "Redux";

    // config strings
    private static final String INPUT_BAM = "input_bam";
    private static final String OUTPUT_BAM = "output_bam";
    private static final String READ_OUTPUTS = "read_output";
    private static final String FORM_CONSENSUS = "form_consensus";
    private static final String READ_LENGTH = "read_length";

    private static final String WRITE_STATS = "write_stats";
    private static final String DROP_DUPLICATES = "drop_duplicates";
    private static final String JITTER_MSI_ONLY = "jitter_msi_only";
    private static final String PARTIION_THREAD_RATIO = "partition_ratio";
    private static final String PARALLEL_CONCATENATION = "parallel_concat";
    private static final String SKIP_FULL_UNMAPPED_READS = "skip_fully_unmapped";
    private static final String SKIP_UNMAPPING = "skip_unmapping";
    private static final String FAIL_SUPP_NO_MATE_CIGAR = "fail_supp_no_mate_cigar";
    private static final String UNMAP_MITOCHONDRIAL = "unmap_mt";

    // dev and options
    public static final String KEEP_INTERIM_BAMS = "keep_interim_bams";
    private static final String NO_WRITE_BAM = "no_write_bam";
    private static final String RUN_CHECKS = "run_checks";
    private static final String SPECIFIC_REGION_FILTER_TYPE = "specific_region_filter";
    private static final String WRITE_READ_BASE_LENGTH = "write_read_base_length";

    public ReduxConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;
        SampleId = configBuilder.getValue(SAMPLE);

        String bamFiles = configBuilder.getValue(INPUT_BAM);
        BamFiles = Arrays.stream(bamFiles.split(CONFIG_FILE_DELIM, -1)).collect(Collectors.toList());

        if(BamFiles.isEmpty())
        {
            RD_LOGGER.error("no BAM files configured");
            System.exit(1);
        }

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = new CachedRefGenome(loadRefGenome(RefGenomeFile));

        SEQUENCING_TYPE = SequencingType.valueOf(configBuilder.getValue(SEQUENCING_TYPE_CFG));

        OutputBam = configBuilder.getValue(OUTPUT_BAM);

        if(OutputBam != null && BamFiles.stream().anyMatch(x -> x.equals(OutputBam)))
        {
            RD_LOGGER.error("output BAM({}) matches input BAM filename", OutputBam);
            System.exit(1);
        }

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
            checkCreateOutputDir(OutputDir);

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

        BQR = new BqrConfig(configBuilder);

        // MD_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        RD_LOGGER.info("output({})", OutputDir);

        BamStringency = BamUtils.validationStringency(configBuilder);

        mReadLength = configBuilder.getInteger(READ_LENGTH);

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        ParallelConcatenation = configBuilder.hasFlag(PARALLEL_CONCATENATION);
        SkipFullyUnmappedReads = configBuilder.hasFlag(SKIP_FULL_UNMAPPED_READS);
        FailOnMissingSuppMateCigar = configBuilder.hasFlag(FAIL_SUPP_NO_MATE_CIGAR);

        DuplicateGroupCollapse = DuplicateGroupCollapseConfig.from(SEQUENCING_TYPE, configBuilder);

        UMIs = UmiConfig.from(configBuilder);

        JitterMsiOnly = configBuilder.hasFlag(JITTER_MSI_ONLY);
        JitterConfig = JitterAnalyserConfig.create(
                SampleId, RefGenomeFile, RefGenVersion, SEQUENCING_TYPE, UMIs.Enabled && UMIs.Duplex, OutputDir, configBuilder);

        FormConsensus = UMIs.Enabled || configBuilder.hasFlag(FORM_CONSENSUS);

        SkipUnmapping = configBuilder.hasFlag(SKIP_UNMAPPING);

        if(configBuilder.hasValue(UNMAP_REGIONS_FILE))
        {
            UnmapRegions = new ReadUnmapper(configBuilder.getValue(UNMAP_REGIONS_FILE));

            if(configBuilder.hasFlag(UNMAP_MITOCHONDRIAL))
                UnmapRegions.addMitochondrialRegion(RefGenVersion);
        }
        else
        {
            Map<String,List<UnmappingRegion>> unmapRegionsMap;

            if(JitterMsiOnly || SkipUnmapping)
            {
                unmapRegionsMap = Collections.emptyMap();
            }
            else
            {
                unmapRegionsMap = Maps.newHashMap();
                ChrBaseRegion excludedRegion = ExcludedRegions.getPolyGRegion(RefGenVersion);
                unmapRegionsMap.put(excludedRegion.Chromosome, Lists.newArrayList(UnmappingRegion.from(excludedRegion, UNMAP_MIN_HIGH_DEPTH)));
            }

            UnmapRegions = new ReadUnmapper(unmapRegionsMap);
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
        PartitionThreadRatio = Threads <= 1 ? 1 : configBuilder.getInteger(PARTIION_THREAD_RATIO);

        LogReadType = ReadOutput.valueOf(configBuilder.getValue(READ_OUTPUTS, NONE.toString()));

        WriteBam = !configBuilder.hasFlag(NO_WRITE_BAM) && !JitterMsiOnly;
        MultiBam = WriteBam && Threads > 1; // now on automatically
        KeepInterimBams = configBuilder.hasFlag(KEEP_INTERIM_BAMS);

        LogReadIds = parseLogReadIds(configBuilder);

        WriteStats = configBuilder.hasFlag(WRITE_STATS);
        PerfDebugTime = configBuilder.getDecimal(PERF_LOG_TIME);
        RunChecks = configBuilder.hasFlag(RUN_CHECKS);
        DropDuplicates = configBuilder.hasFlag(DROP_DUPLICATES);
        WriteReadBaseLength = configBuilder.getInteger(WRITE_READ_BASE_LENGTH);

        if(RunChecks)
        {
            RD_LOGGER.info("running debug options: read-checks({})", RunChecks);
        }

        mReadChecker = new ReadChecker(RunChecks);
    }

    // convenience
    public static boolean isIllumina() { return SEQUENCING_TYPE == ILLUMINA; }
    public static boolean isSbx() { return SEQUENCING_TYPE == SBX; }
    public static boolean isUltima() { return SEQUENCING_TYPE == ULTIMA; }
    public static boolean isBiomdal() { return SEQUENCING_TYPE == BIOMODAL; }

    public boolean isValid() { return mIsValid; }

    public int readLength() { return mReadLength; }
    public void setReadLength(int readLength)
    {
        mReadLength = readLength;
    }

    public boolean perfDebug() { return PerfDebugTime > 0; }
    public ReadChecker readChecker() { return mReadChecker; }

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
        configBuilder.addPaths(INPUT_BAM, false, "BAM file path, separated by ',' if multiple");
        configBuilder.addConfigItem(OUTPUT_BAM, false, "Output BAM filename");
        addRefGenomeConfig(configBuilder, true);
        SequencingType.registerConfig(configBuilder);
        configBuilder.addInteger(READ_LENGTH, "Read length, otherwise will sample from BAM", 0);

        BqrConfig.registerConfig(configBuilder);

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

        JitterAnalyserConfig.addConfig(configBuilder);
        DuplicateGroupCollapseConfig.addConfig(configBuilder);

        addThreadOptions(configBuilder);
        configBuilder.addInteger(PARTIION_THREAD_RATIO, "Partitions per thread, impacts BAM-writing performance", 2);
        configBuilder.addFlag(PARALLEL_CONCATENATION, "Parallel final BAM concatenation");
        configBuilder.addFlag(SKIP_FULL_UNMAPPED_READS, "Skip processing existing fully unmapped reads");
        configBuilder.addFlag(SKIP_UNMAPPING, "Skip unmapping routine, including excluded regions");
        configBuilder.addFlag(UNMAP_MITOCHONDRIAL, "Unmap mitochondrial reads");

        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, PERF_LOG_TIME_DESC, 0);
        configBuilder.addFlag(RUN_CHECKS, "Run duplicate mismatch checks");
        configBuilder.addFlag(FAIL_SUPP_NO_MATE_CIGAR, "Fail if supplementary is missing mate CIGAR ");
        configBuilder.addConfigItem(SPECIFIC_REGION_FILTER_TYPE, "Used with specific regions, to filter mates or supps");

        configBuilder.addInteger(WRITE_READ_BASE_LENGTH, "Number of read bases to write with read data", 0);
    }

    @VisibleForTesting
    public ReduxConfig(final RefGenomeInterface refGenome, boolean umiEnabled, boolean duplexUmi, boolean formConsensus,
            final ReadUnmapper readUnmapper, final SequencingType sequencingType, int sbxMaxDuplicateDistance)
    {
        mIsValid = true;
        SampleId = "";
        BamFiles = Lists.newArrayList();
        RefGenomeFile = null;
        OutputBam = null;
        OutputDir = null;
        OutputId = "";
        RefGenVersion = V37;
        RefGenome = new CachedRefGenome(refGenome);

        BamStringency = ValidationStringency.STRICT;
        mReadLength = DEFAULT_READ_LENGTH;

        UMIs = new UmiConfig(umiEnabled, duplexUmi, String.valueOf(DEFAULT_DUPLEX_UMI_DELIM), false);
        FormConsensus = umiEnabled || formConsensus;

        BQR = new BqrConfig();

        SpecificChrRegions = new SpecificRegions();
        SpecificRegionsFilterType = FilterReadsType.MATE_AND_SUPP;

        BamToolPath = null;
        ParallelConcatenation = false;

        UnmapRegions = readUnmapper;

        JitterMsiOnly = false;
        JitterConfig = null;

        DuplicateGroupCollapse = new DuplicateGroupCollapseConfig(sequencingType, sbxMaxDuplicateDistance);

        WriteBam = false;
        MultiBam = false;
        KeepInterimBams = false;
        SkipFullyUnmappedReads = false;
        SkipUnmapping = false;
        LogReadType = NONE;
        FailOnMissingSuppMateCigar = false;

        LogReadIds = Lists.newArrayList();
        Threads = 0;
        PartitionThreadRatio = 1;
        PerfDebugTime = 0;
        RunChecks = true;
        WriteStats = false;
        DropDuplicates = false;
        WriteReadBaseLength = 0;

        mReadChecker = new ReadChecker(false);
    }

    @VisibleForTesting
    public ReduxConfig(final RefGenomeInterface refGenome, boolean umiEnabled, boolean duplexUmi, boolean formConsensus,
            final ReadUnmapper readUnmapper, final SequencingType sequencingType)
    {
        this(refGenome, umiEnabled, duplexUmi, formConsensus, readUnmapper, sequencingType, 0);
    }

    @VisibleForTesting
    public ReduxConfig(
            final RefGenomeInterface refGenome, boolean umiEnabled, boolean duplexUmi, boolean formConsensus,
            final ReadUnmapper readUnmapper)
    {
        this(refGenome, umiEnabled, duplexUmi, formConsensus, readUnmapper, ILLUMINA);
    }
}
