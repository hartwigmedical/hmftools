package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE_DESC;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.REGIONS_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.loadSpecificRegionsConfig;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;

import java.util.List;

import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class CheckConfig
{
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;

    public final String OutputDir;
    public final String OutputId;
    public final String OutputBam;

    public final int Threads;
    public final String BamToolPath;
    public final boolean SkipUnmapped;
    public final boolean WriteIncompleteFragments;
    public final boolean DropIncompleteFragments;

    // debug
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final boolean PerfDebug;

    public static final String OUTPUT_BAM_FILE = "output_bam";
    public static final String SKIP_UNMAPPED = "skip_unmapped";
    public static final String WRITE_INCOMPLETE_FRAGS = "write_incompletes";
    public static final String DROP_INCOMPLETE_FRAGS = "drop_incompletes";
    public static final String CFG_LOG_READ_COUNT = "log_read_count";

    protected static int LOG_READ_COUNT = 10_000_000;

    public CheckConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        if(configBuilder.hasValue(OUTPUT_BAM_FILE))
        {
            OutputBam = configBuilder.getValue(OUTPUT_BAM_FILE);
            OutputDir = pathFromFile(OutputBam);
        }
        else
        {
            OutputDir = configBuilder.hasValue(OUTPUT_DIR) ? parseOutputDir(configBuilder) : pathFromFile(BamFile);
            OutputBam = null;
        }

        OutputId = configBuilder.getValue(OUTPUT_ID);

        RefGenVersion = deriveRefGenomeVersion(BamFile);
        BT_LOGGER.info("input bam({}) refGenome({})", BamFile, RefGenVersion);
        BT_LOGGER.info("output({})", OutputDir);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        SpecificChrRegions = new SpecificRegions();
        loadSpecificRegionsConfig(configBuilder, SpecificChrRegions.Chromosomes, SpecificChrRegions.Regions);

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        Threads = parseThreads(configBuilder);
        LogReadIds = parseLogReadIds(configBuilder);

        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
        SkipUnmapped = configBuilder.hasFlag(SKIP_UNMAPPED);
        WriteIncompleteFragments = configBuilder.hasFlag(WRITE_INCOMPLETE_FRAGS);
        DropIncompleteFragments = configBuilder.hasFlag(DROP_INCOMPLETE_FRAGS);

        LOG_READ_COUNT = configBuilder.getInteger(CFG_LOG_READ_COUNT);
    }

    public String formFilename(final String fileId, final String fileExtension)
    {
        String bamFilename = filenamePart(BamFile);
        String bamPrefix = bamFilename.substring(0, bamFilename.lastIndexOf("."));
        String outputFile = OutputDir + bamPrefix + "." + fileId;

        if(OutputId != null)
        {
            outputFile += "." + OutputId;
        }

        outputFile += fileExtension;
        return outputFile;
    }

    public boolean writeBam() { return BamToolPath != null; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, BAM_FILE_DESC);

        configBuilder.addPath(
                REGIONS_FILE, false,
                "TSV or BED file with regions to slice, expected columns Chromosome,PositionStart,PositionEnd or no headers");

        addRefGenomeFile(configBuilder, true);;

        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addInteger(CFG_LOG_READ_COUNT, "Log partition processed read count frequency", LOG_READ_COUNT);
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addFlag(SKIP_UNMAPPED, "Skip full unmapped reads");
        configBuilder.addFlag(WRITE_INCOMPLETE_FRAGS, "Write incomplete fragments to TSV");
        configBuilder.addFlag(DROP_INCOMPLETE_FRAGS, "Drop incomplete fragments from BAM");

        configBuilder.addConfigItem(OUTPUT_BAM_FILE, "Output BAM filename, default is input BAM + 'final'");

        BamToolName.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
    }
}
