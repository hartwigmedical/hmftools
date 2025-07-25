package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE_DESC;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.REGIONS_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.checkFileExists;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.loadSpecificRegionsConfig;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegionList;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;

import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SliceConfig
{
    public final String BamFile;
    public final String OutputPrefix;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;

    public final String OutputDir;

    public final boolean WriteBam;
    public final boolean UnsortedBam;
    public final boolean WriteReads;
    public final boolean DropExcluded;
    public final boolean OnlySupplementaries;
    public final boolean WriteReadBases;
    public final boolean LogMissingReads;
    public final int MaxRemoteReads;
    public final int MaxPartitionReads;
    public final int MaxUnmappedReads;
    public final int Threads;
    public final String BamToolPath;

    // debug
    public final SpecificRegions SliceRegions;
    public final List<String> LogReadIds;
    public final boolean PerfDebug;

    private boolean mIsValid;

    private static final String OUTPUT_PREFIX = "output_prefix";
    private static final String WRITE_BAM = "write_bam";
    private static final String UNSORTED_BAM = "unsorted_bam";
    private static final String WRITE_READS = "write_reads";
    private static final String DROP_EXCLUDED = "drop_excluded";
    private static final String ONLY_SUPPS = "only_supps";
    private static final String MAX_PARTITION_READS = "max_partition_reads";
    private static final String MAX_REMOTE_READS = "max_remote_reads";
    private static final String WRITE_READ_BASES = "write_read_bases";
    private static final String MAX_UNMAPPED_READS = "max_unmapped_reads";
    private static final String LOG_MISSING_READS = "log_missing_reads";

    public static final int UNMAPPED_READS_DISABLED = -1;
    public static final int UNMAPPED_READS_ALL = 0;

    public SliceConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        BamFile = configBuilder.getValue(BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        OutputDir = configBuilder.hasValue(OUTPUT_DIR) ? parseOutputDir(configBuilder) : pathFromFile(BamFile);

        if(configBuilder.hasValue(OUTPUT_PREFIX))
        {
            OutputPrefix = configBuilder.getValue(OUTPUT_PREFIX);
        }
        else
        {
            String filename = filenamePart(BamFile);
            int extIndex = filename.lastIndexOf(".");
            String bamFileName = filename.substring(0, extIndex);

            OutputPrefix = bamFileName + ".slice";
        }

        WriteReads = configBuilder.hasFlag(WRITE_READS);
        WriteBam = !WriteReads || configBuilder.hasFlag(WRITE_BAM) || configBuilder.hasFlag(BAMTOOL_PATH);
        UnsortedBam = configBuilder.hasFlag(UNSORTED_BAM);
        DropExcluded = configBuilder.hasFlag(DROP_EXCLUDED);
        WriteReadBases = WriteReads && configBuilder.hasFlag(WRITE_READ_BASES);
        OnlySupplementaries = configBuilder.hasFlag(ONLY_SUPPS);
        LogMissingReads = configBuilder.hasFlag(LOG_MISSING_READS);
        MaxRemoteReads = configBuilder.getInteger(MAX_REMOTE_READS);
        MaxPartitionReads = configBuilder.getInteger(MAX_PARTITION_READS);
        MaxUnmappedReads = configBuilder.getInteger(MAX_UNMAPPED_READS);

        if(BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: bam({}) refGenome({}) outputDir({})",
                    BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = deriveRefGenomeVersion(BamFile);
        BT_LOGGER.info("input bam({}) refGenome({})", BamFile, RefGenVersion);
        BT_LOGGER.info("output({}) and file prefix({})", OutputDir, OutputPrefix);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        SliceRegions = new SpecificRegions();

        if(configBuilder.hasValue(REGIONS_FILE))
        {
            List<ChrBaseRegion> regions = loadChrBaseRegionList(configBuilder.getValue(REGIONS_FILE));
            regions.forEach(x -> SliceRegions.addRegion(x));
        }
        else
        {
            mIsValid &= loadSpecificRegionsConfig(configBuilder, SliceRegions.Chromosomes, SliceRegions.Regions);
        }

        if(SliceRegions.Regions.isEmpty())
        {
            BT_LOGGER.error("missing specific regions or slice BED file for slicing");
            System.exit(1);
        }

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        Threads = parseThreads(configBuilder);

        LogReadIds = parseLogReadIds(configBuilder);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        mIsValid = checkFileExists(BamFile) && checkFileExists(RefGenomeFile);
        return mIsValid;
    }

    public String formFilename(final WriteType fileType)
    {
        String outputFile = OutputDir + OutputPrefix + fileType.extension();
        return outputFile;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, BAM_FILE_DESC);

        configBuilder.addPath(
                REGIONS_FILE, false,
                "TSV or BED file with regions to slice, expected columns Chromosome,PositionStart,PositionEnd or no headers");

        addRefGenomeFile(configBuilder, true);;
        configBuilder.addConfigItem(OUTPUT_PREFIX, "File prefix for BAM and read TSV");

        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addInteger(MAX_PARTITION_READS, "Max partition reads (perf-only)", 0);
        configBuilder.addInteger(MAX_REMOTE_READS, "Max remote reads (perf-only)", 0);
        configBuilder.addInteger(MAX_UNMAPPED_READS, "Max unmapped reads: 0 means all, -1 to disable (default)", -1);
        configBuilder.addFlag(WRITE_BAM, "Write BAM file for sliced region");
        configBuilder.addFlag(UNSORTED_BAM, "Write BAM unsorted");
        configBuilder.addFlag(WRITE_READS, "Write reads file for sliced region");
        configBuilder.addFlag(WRITE_READ_BASES, "Write read bases to TSV file");
        configBuilder.addFlag(DROP_EXCLUDED, "Ignore remote reads in excluded regions (eg poly-G)");
        configBuilder.addFlag(ONLY_SUPPS, "Only capture supplementary reads");
        configBuilder.addFlag(LOG_MISSING_READS, "Log missing reads");
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        BamToolName.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
    }

    @VisibleForTesting
    public SliceConfig()
    {
        mIsValid = true;
        OutputPrefix = "";
        BamFile = "";
        RefGenomeFile = "";
        RefGenVersion = V37;
        OutputDir = "";
        WriteReads = false;
        WriteBam = false;
        WriteReadBases = false;
        UnsortedBam = false;
        DropExcluded = false;
        OnlySupplementaries = false;
        LogMissingReads = false;
        MaxRemoteReads = 0;
        MaxPartitionReads = 0;
        MaxUnmappedReads = UNMAPPED_READS_DISABLED;
        PartitionSize = DEFAULT_CHR_PARTITION_SIZE;
        SliceRegions = new SpecificRegions();
        Threads = 0;
        BamToolPath = null;
        PerfDebug = false;
        LogReadIds = Collections.emptyList();
    }
}
