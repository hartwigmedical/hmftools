package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.Level;

public class CompareConfig
{
    public final String OutputFile;
    public final String OrigBamFile;
    public final String NewBamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;
    public final int MaxCachedReadsPerThread;
    public final boolean ExcludeRegions;
    public final boolean IgnoreDupDiffs;
    public final boolean IgnoreAlterations; // consensus reads and internal unmappings

    public final boolean IgnoreConsensusReads;

    public final boolean IgnoreSupplementaryReads;

    public final boolean MatchOrigUnmapped;
    public final boolean MatchNewUnmapped;

    public final int Threads;
    public final List<String> LogReadIds;

    // debug
    public final SpecificRegions SpecificChrRegions;

    private static final String OUTPUT_FILE = "output_file";
    private static final String ORIG_BAM_FILE = "orig_bam_file";
    private static final String NEW_BAM_FILE = "new_bam_file";
    private static final String EXCLUDE_REGIONS = "exclude_regions";
    private static final String MAX_CACHED_READS_PER_THREAD = "max_cached_reads_per_thread";
    private static final String IGNORE_DUP_DIFFS = "ignore_dup_diffs";
    private static final String IGNORE_ALTERATIONS = "ignore_alterations";
    private static final String IGNORE_SUPPLEMENTARY_READS = "ignore_supplementary_reads";
    private static final String IGNORE_CONSENSUS_READS = "ignore_consensus_reads";
    private static final String MATCH_ORIG_UNMAPPED = "match_orig_unmapped";
    private static final String MATCH_NEW_UNMAPPED = "match_new_unmapped";

    private static final int DEFAULT_CHR_PARTITION_SIZE = 10_000_000;

    public CompareConfig(final ConfigBuilder configBuilder)
    {
        OutputFile =  configBuilder.getValue(OUTPUT_FILE);
        OrigBamFile =  configBuilder.getValue(ORIG_BAM_FILE);
        NewBamFile =  configBuilder.getValue(NEW_BAM_FILE);
        RefGenomeFile =  configBuilder.getValue(REF_GENOME);

        if(OrigBamFile == null || NewBamFile == null || OutputFile == null)
        {
            BT_LOGGER.error("missing config: bam(orig={} new={}) outputDir({})",
                    OrigBamFile != null, NewBamFile != null, OutputFile != null);
            System.exit(1);
        }

        RefGenVersion = BamUtils.deriveRefGenomeVersion(OrigBamFile);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        int maxCachedReadsPerThread = configBuilder.getInteger(MAX_CACHED_READS_PER_THREAD);
        IgnoreDupDiffs = configBuilder.hasFlag(IGNORE_DUP_DIFFS);
        IgnoreAlterations = configBuilder.hasFlag(IGNORE_ALTERATIONS);
        IgnoreConsensusReads = configBuilder.hasFlag(IGNORE_CONSENSUS_READS);
        IgnoreSupplementaryReads = configBuilder.hasFlag(IGNORE_SUPPLEMENTARY_READS);

        MatchOrigUnmapped = configBuilder.hasFlag(MATCH_ORIG_UNMAPPED);
        MatchNewUnmapped = configBuilder.hasFlag(MATCH_NEW_UNMAPPED);

        BT_LOGGER.info("refGenomeVersion({}) origBam({}) newBam({})", RefGenVersion, OrigBamFile, NewBamFile);
        BT_LOGGER.info("origBam({}) newBam({})", OrigBamFile, NewBamFile);
        BT_LOGGER.info("outputFile({})", OutputFile);
        BT_LOGGER.info("ignoreDupDiffs({}) ignoreAlterations({}) ignoreConsensusReads({}) ignoreSupplementaryReads({})",
                IgnoreDupDiffs, IgnoreAlterations, IgnoreConsensusReads, IgnoreSupplementaryReads);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            System.exit(1);

        ExcludeRegions = configBuilder.hasFlag(EXCLUDE_REGIONS);

        Threads = Math.max(parseThreads(configBuilder), 1);

        LogReadIds = parseLogReadIds(configBuilder);

        if(maxCachedReadsPerThread == 0)
        {
            // calculate appropriate limit based on how much memory is given
            int maxMemKB = (int)(Runtime.getRuntime().maxMemory() / 1024);
            // seems 5kB for 1 read is about right
            MaxCachedReadsPerThread = maxMemKB / 5 / Threads;
        }
        else
        {
            MaxCachedReadsPerThread = maxCachedReadsPerThread;
        }

        BT_LOGGER.printf(Level.INFO, "maxCachedReadsPerThread(%,d)", MaxCachedReadsPerThread);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output comparison file");
        configBuilder.addPath(ORIG_BAM_FILE, true, "Original BAM file");
        configBuilder.addPath(NEW_BAM_FILE,true, "New BAM file");
        configBuilder.addInteger(PARTITION_SIZE, "Partition size (performance only)", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addInteger(MAX_CACHED_READS_PER_THREAD, "Maximum cached reads per thread", 0);

        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addFlag(EXCLUDE_REGIONS, "If set, ignore excluded regions");
        configBuilder.addFlag(IGNORE_DUP_DIFFS, "If set, ignore duplicate diffs");
        configBuilder.addFlag(IGNORE_ALTERATIONS, "If set, ignore consensus reads and internal unmappings");
        configBuilder.addFlag(IGNORE_CONSENSUS_READS, "If set, ignore consensus reads");
        configBuilder.addFlag(IGNORE_SUPPLEMENTARY_READS, "If set, ignore supplementary reads");

        configBuilder.addFlag(MATCH_ORIG_UNMAPPED, "If set, match unmapped reads in original bam against corresponding read in new bam");
        configBuilder.addFlag(MATCH_NEW_UNMAPPED, "If set, match unmapped reads in new bam against corresponding read in original bam");

        addRefGenomeFile(configBuilder, false);
        addSpecificChromosomesRegionsConfig(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    public boolean ignoreUnmapped() { return SpecificChrRegions != null && !SpecificChrRegions.Regions.isEmpty(); }

    @VisibleForTesting
    public CompareConfig()
    {
        OutputFile = null;
        OrigBamFile = null;
        NewBamFile = null;
        RefGenomeFile = null;
        RefGenVersion = null;
        PartitionSize = 0;
        MaxCachedReadsPerThread = 0;
        ExcludeRegions = false;
        IgnoreDupDiffs = false;
        IgnoreAlterations = false;
        IgnoreConsensusReads = false;
        IgnoreSupplementaryReads = false;
        MatchOrigUnmapped = false;
        MatchNewUnmapped = false;
        Threads = 0;
        LogReadIds = null;
        SpecificChrRegions = null;
    }
}
