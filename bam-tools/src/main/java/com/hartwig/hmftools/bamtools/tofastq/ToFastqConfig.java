package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE_DESC;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public class ToFastqConfig
{
    public final String BamFile;

    @Nullable public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final String OutputDir;
    @Nullable public final String OutputId;
    public final FileSplitMode SplitMode;
    public final int Threads;
    public final int PartitionSize;

    public final SpecificRegions SpecificChrRegions;
    public final boolean PerfDebug;

    private static final String FILE_SPLIT_MODE = "split_mode";

    public static final String CHR_UNMAPPED = "unmapped"; // to test unmapped reads

    private static final int DEFAULT_PARTITION_SIZE = 1_000_000;

    public ToFastqConfig(final ConfigBuilder configBuilder)
    {
        BamFile =  configBuilder.getValue(BAM_FILE);
        RefGenomeFile =  configBuilder.getValue(REF_GENOME);

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
            checkCreateOutputDir(OutputDir);
        }
        else
        {
            OutputDir = pathFromFile(BamFile);
        }

        OutputId =  configBuilder.getValue(OUTPUT_ID);

        if(BamFile == null || OutputDir == null)
        {
            BT_LOGGER.error("missing config: bam({}) outputDir({})",
                    BamFile != null, OutputDir != null);
            System.exit(1);
        }

        RefGenVersion = deriveRefGenomeVersion(BamFile);
        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        BT_LOGGER.info("output({})", OutputDir);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        SpecificChrRegions = SpecificRegions.from(configBuilder);
        if(SpecificChrRegions == null)
        {
            System.exit(1);
        }

        Threads = Math.max(parseThreads(configBuilder), 1);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
        SplitMode = FileSplitMode.valueOf(configBuilder.getValue(FILE_SPLIT_MODE));

        BT_LOGGER.info("threads({})", Threads);
        BT_LOGGER.info("splitMode({})", SplitMode);
    }

    public String formFilePrefix(final String threadId, final String readGroupId, boolean includeBamFilePrefix)
    {
        List<String> filenameParts = new ArrayList<>();

        if(includeBamFilePrefix)
        {
            String bamFile = filenamePart(BamFile);
            filenameParts.add(bamFile.substring(0, bamFile.indexOf(".")));
        }

        if(!threadId.isEmpty())
            filenameParts.add(threadId);

        if(!readGroupId.isEmpty())
            filenameParts.add(readGroupId);

        if(OutputId != null)
            filenameParts.add(OutputId);

        return OutputDir + String.join(".", filenameParts);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, false);

        configBuilder.addPath(BAM_FILE, true, BAM_FILE_DESC);

        configBuilder.addConfigItem(FILE_SPLIT_MODE, false, "File split mode, NONE, READ_GROUP (default), THREAD",
                FileSplitMode.READ_GROUP.name());
        configBuilder.addInteger(PARTITION_SIZE, "Partition split size", DEFAULT_PARTITION_SIZE);
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
