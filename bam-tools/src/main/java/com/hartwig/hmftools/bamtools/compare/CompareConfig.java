package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class CompareConfig
{
    public final String OutputFile;
    public final String RefBamFile;
    public final String NewBamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;
    public final int MaxPartitionReads;
    public final boolean ExcludeRegions;
    public final boolean IgnoreDupDiffs;
    public final boolean IgnoreAlterations; // consensus reads and internal unmappings

    public final int Threads;
    public final List<String> LogReadIds;

    // debug
    public final SpecificRegions SpecificChrRegions;

    private static final String OUTPUT_FILE = "output_file";
    private static final String REF_BAM_FILE = "ref_bam_file";
    private static final String NEW_BAM_FILE = "new_bam_file";
    private static final String EXCLUDE_REGIONS = "exclude_regions";
    private static final String MAX_PARTITION_READS = "max_partition_reads";
    private static final String IGNORE_DUP_DIFFS = "ignore_dup_diffs";
    private static final String IGNORE_ALTERATIONS = "ignore_alterations";

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

        BT_LOGGER.info("refBam({}) newBam({})", RefBamFile, NewBamFile);
        BT_LOGGER.info("output file({})", OutputFile);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        MaxPartitionReads = configBuilder.getInteger(MAX_PARTITION_READS);
        IgnoreDupDiffs = configBuilder.hasFlag(IGNORE_DUP_DIFFS);
        IgnoreAlterations = configBuilder.hasFlag(IGNORE_ALTERATIONS);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            System.exit(1);

        ExcludeRegions = configBuilder.hasFlag(EXCLUDE_REGIONS);

        Threads = parseThreads(configBuilder);

        LogReadIds = parseLogReadIds(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addInteger(MAX_PARTITION_READS, "Maximum partition reads before exit", 0);

        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output comparison file");
        configBuilder.addRequiredConfigItem(REF_BAM_FILE, "Ref BAM file");
        configBuilder.addRequiredConfigItem(NEW_BAM_FILE,"New BAM file");
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addFlag(EXCLUDE_REGIONS, "Ignore excluded regions");
        configBuilder.addFlag(IGNORE_DUP_DIFFS, "Ignore duplicate diffs");
        configBuilder.addFlag(IGNORE_ALTERATIONS, "Ignore consensus reads and internal unmappings");

        addRefGenomeConfig(configBuilder, true);;
        addSpecificChromosomesRegionsConfig(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
