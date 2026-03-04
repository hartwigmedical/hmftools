package com.hartwig.hmftools.bamtools.utils;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
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
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class UmiAnalyserConfig
{
    public final String UmiBamFile;
    public final String RawBamFile;
    public final String KnownUmiFile;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final String OutputDir;
    public final String OutputId;

    public final int PartitionSize;
    public final int Threads;

    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;

    private static final String UMI_BAM_FILE = "umi_bam_file";
    private static final String RAW_BAM_FILE = "raw_bam_file";
    private static final String KNOWN_UMI_FILE = "known_umi_file";

    public UmiAnalyserConfig(final ConfigBuilder configBuilder)
    {
        UmiBamFile = configBuilder.getValue(UMI_BAM_FILE);
        RawBamFile = configBuilder.getValue(RAW_BAM_FILE);
        KnownUmiFile = configBuilder.getValue(KNOWN_UMI_FILE);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        RefGenVersion = BamUtils.deriveRefGenomeVersion(RawBamFile);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        Threads = parseThreads(configBuilder);

        LogReadIds = parseLogReadIds(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(UMI_BAM_FILE, true, "BAM with UMIs extracted");
        configBuilder.addConfigItem(RAW_BAM_FILE, true, "Raw BAM without UMIs extracted");
        configBuilder.addPath(KNOWN_UMI_FILE, false, "File with known UMI sequences");

        addRefGenomeFile(configBuilder, true);
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);

        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, PERF_LOG_TIME_DESC, 0);
    }
}
