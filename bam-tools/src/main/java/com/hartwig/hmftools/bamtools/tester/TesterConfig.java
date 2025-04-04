package com.hartwig.hmftools.bamtools.tester;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE_DESC;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SEQUENCING_TYPE_CFG;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_LOG_TIME;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_LOG_TIME_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import htsjdk.samtools.ValidationStringency;

public class TesterConfig
{
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    // public final RefGenomeInterface RefGenome;
    public final SequencingType Sequencing;

    public final ValidationStringency BamStringency;
    public final String OutputDir;
    public final String OutputId;

    public final int PartitionSize;
    public final int Threads;

    public final String BamToolPath;
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final double PerfDebugTime;


    public TesterConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        Sequencing = SequencingType.parseConfig(configBuilder);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        RefGenVersion = BamUtils.deriveRefGenomeVersion(BamFile);

        BamStringency = BamUtils.validationStringency(configBuilder);
        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        Threads = parseThreads(configBuilder);

        LogReadIds = parseLogReadIds(configBuilder);

        PerfDebugTime = configBuilder.getDecimal(PERF_LOG_TIME);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(BAM_FILE, true, BAM_FILE_DESC);
        addRefGenomeFile(configBuilder, true);
        SequencingType.registerConfig(configBuilder);
        BamToolName.addConfig(configBuilder);
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);

        addValidationStringencyOption(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, PERF_LOG_TIME_DESC, 0);
    }
}
