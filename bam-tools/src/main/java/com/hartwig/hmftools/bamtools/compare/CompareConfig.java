package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.addCommonCommandOptions;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.loadSpecificRegionsConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class CompareConfig
{
    public final String OutputFile;
    public final String RefBamFile;
    public final String NewBamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;

    public final int Threads;
    public final List<String> LogReadIds;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<ChrBaseRegion> SpecificRegions;

    private static final String OUTPUT_FILE = "output_file";
    private static final String REF_BAM_FILE = "ref_bam_file";
    private static final String NEW_BAM_FILE = "new_bam_file";

    private static final int DEFAULT_CHR_PARTITION_SIZE = 100000;

    public CompareConfig(final ConfigBuilder configBuilder)
    {
        OutputFile =  configBuilder.getValue(OUTPUT_FILE);
        RefBamFile =  configBuilder.getValue(REF_BAM_FILE);
        NewBamFile =  configBuilder.getValue(NEW_BAM_FILE);
        RefGenomeFile =  configBuilder.getValue(REF_GENOME);

        if(RefBamFile == null || NewBamFile == null || OutputFile == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: bam(ref={} new={}) refGenome({}) outputDir({})",
                    RefBamFile != null, NewBamFile != null, RefGenomeFile != null, OutputFile != null);
            System.exit(1);
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, RefBamFile);
        BT_LOGGER.info("output file({})", OutputFile);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        if(!loadSpecificRegionsConfig(configBuilder, SpecificChromosomes, SpecificRegions))
            System.exit(1);

        Threads = parseThreads(configBuilder);

        LogReadIds = parseLogReadIds(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addCommonCommandOptions(configBuilder);

        configBuilder.addIntegerItem(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);

        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output comparison file");
        configBuilder.addRequiredConfigItem(REF_BAM_FILE, "Ref BAM file");
        configBuilder.addRequiredConfigItem(NEW_BAM_FILE,"New BAM file");
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
    }
}
