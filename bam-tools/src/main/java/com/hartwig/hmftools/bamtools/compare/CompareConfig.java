package com.hartwig.hmftools.bamtools.compare;

import static java.lang.Math.max;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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
    public final boolean IgnoreConsensusReads;
    public final boolean IgnoreSupplementaryReads;
    public final boolean IgnoreSecondaryReads;
    public final boolean IgnoreSupplementaryAttribute;
    public final boolean IgnoreReduxUnmapped;
    public final boolean CompareCoordsOnly;
    public final boolean CheckBasesAndQuals;
    public final boolean StandardChromosomes;
    public final int CigarBoundaryTolerance;

    public final boolean IgnoreReduxAlterations; // consensus reads and internal unmappings

    public final int Threads;
    public final List<String> LogReadIds;

    // comma-separated list of TSV paths (Chromosome / PosStart / PosEnd headers; extra columns ignored).
    // any readname whose any alignment overlaps an ignore region by > 50% of its reference span is dropped
    // from comparison in both BAMs.
    public final List<String> IgnoreRegionFiles;

    // debug
    public final SpecificRegions SpecificChrRegions;

    private static final String OUTPUT_FILE = "output_file";
    private static final String ORIG_BAM_FILE = "orig_bam_file";
    private static final String NEW_BAM_FILE = "new_bam_file";
    private static final String EXCLUDE_REGIONS = "exclude_regions";
    private static final String MAX_CACHED_READS_PER_THREAD = "max_cached_reads_per_thread";

    private static final String IGNORE_DUP_DIFFS = "ignore_dup_diffs";
    private static final String IGNORE_ALTERATIONS = "ignore_alterations";
    private static final String IGNORE_SUPPLEMENTARY_READS = "ignore_supp_reads";
    private static final String IGNORE_SECONDARY_READS = "ignore_sec_reads";
    private static final String IGNORE_SUPPLEMENTARY_ATTRIBUTE = "ignore_supp_attribute";
    private static final String IGNORE_CONSENSUS_READS = "ignore_consensus_reads";
    private static final String IGNORE_REDUX_UNMAPPED = "ignore_redux_unmapped";
    private static final String IGNORE_REDUX_DIFFS = "ignore_redux_diffs";
    private static final String IGNORE_REGIONS_FILES = "ignore_regions";
    private static final String CIGAR_BOUNDARY_TOLERANCE = "cigar_boundary_tolerance";
    private static final String COORDS_ONLY = "coords_only";
    private static final String CHECK_BASES_QUALS = "check_bases_quals";
    protected static final String STANDARD_CHROMOSOMES = "std_chromosomes";

    private static final int DEFAULT_CHR_PARTITION_SIZE = 10_000_000;

    protected static final String FULLY_UNMAPPED_PARTITION = "fully_unmapped";

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

        CompareCoordsOnly = configBuilder.hasFlag(COORDS_ONLY);

        CigarBoundaryTolerance = configBuilder.getInteger(CIGAR_BOUNDARY_TOLERANCE);

        if(configBuilder.hasFlag(IGNORE_REDUX_DIFFS))
        {
            IgnoreDupDiffs = true;
            IgnoreConsensusReads = true;
            IgnoreSupplementaryReads = true;
            IgnoreSecondaryReads = true;
            IgnoreSupplementaryAttribute = true;
            IgnoreReduxUnmapped = true;
            CheckBasesAndQuals = false;

            // setting these false allows primary read differences other than those caused by Redux to be evaluated - ie they are not dropped
            IgnoreReduxAlterations = false;
        }
        else
        {
            IgnoreDupDiffs = configBuilder.hasFlag(IGNORE_DUP_DIFFS);
            IgnoreConsensusReads = configBuilder.hasFlag(IGNORE_CONSENSUS_READS);
            IgnoreSupplementaryReads = configBuilder.hasFlag(IGNORE_SUPPLEMENTARY_READS);
            IgnoreSecondaryReads = configBuilder.hasFlag(IGNORE_SECONDARY_READS);
            IgnoreSupplementaryAttribute = configBuilder.hasFlag(IGNORE_SUPPLEMENTARY_ATTRIBUTE);
            IgnoreReduxUnmapped = configBuilder.hasFlag(IGNORE_REDUX_UNMAPPED);
            IgnoreReduxAlterations = configBuilder.hasFlag(IGNORE_ALTERATIONS);
            CheckBasesAndQuals = configBuilder.hasFlag(CHECK_BASES_QUALS);
        }

        BT_LOGGER.info("origBam({}) newBam({})", OrigBamFile, NewBamFile);

        BT_LOGGER.info("outputFile({})", OutputFile);

        BT_LOGGER.info("ignore options: duplicateDiffs({}) supps(reads={} attributes={}) Redux(all={} unmapped={} consensus={} basesQuals={})",
                IgnoreDupDiffs, IgnoreSupplementaryReads, IgnoreSupplementaryAttribute,
                IgnoreReduxAlterations, IgnoreReduxUnmapped, IgnoreConsensusReads, CheckBasesAndQuals);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            System.exit(1);

        StandardChromosomes = configBuilder.hasFlag(STANDARD_CHROMOSOMES);

        if(StandardChromosomes)
        {
            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                SpecificChrRegions.Chromosomes.add(RefGenVersion.versionedChromosome(chromosome.toString()));
            }
        }

        ExcludeRegions = configBuilder.hasFlag(EXCLUDE_REGIONS);

        final String ignoreRegionFiles = configBuilder.getValue(IGNORE_REGIONS_FILES);
        if(ignoreRegionFiles != null && !ignoreRegionFiles.isEmpty())
        {
            IgnoreRegionFiles = Arrays.asList(ignoreRegionFiles.split(","));
            BT_LOGGER.info("ignore-region TSV file(s): {}", IgnoreRegionFiles);
        }
        else
        {
            IgnoreRegionFiles = Collections.emptyList();
        }

        Threads = max(parseThreads(configBuilder), 1);

        LogReadIds = parseLogReadIds(configBuilder);

        if(maxCachedReadsPerThread == 0)
        {
            // calculate appropriate limit based on how much memory is given
            int maxMemKB = (int)(Runtime.getRuntime().maxMemory() / 1024);
            // seems 5kB for 1 read is about right
            MaxCachedReadsPerThread = maxMemKB / 5 / Threads;

            BT_LOGGER.debug("setting maxCachedReadsPerThread({})", MaxCachedReadsPerThread);
        }
        else
        {
            MaxCachedReadsPerThread = max(maxCachedReadsPerThread, 0);
        }
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output comparison file");
        configBuilder.addPath(ORIG_BAM_FILE, true, "Original BAM file");
        configBuilder.addPath(NEW_BAM_FILE,true, "New BAM file");
        configBuilder.addInteger(PARTITION_SIZE, "Partition size (performance only)", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addInteger(MAX_CACHED_READS_PER_THREAD, "Maximum cached reads per thread", 0);

        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addFlag(EXCLUDE_REGIONS, "Specify regions to exclude");

        configBuilder.addFlag(COORDS_ONLY, "Only compare read coordinates");
        configBuilder.addFlag(IGNORE_REDUX_DIFFS, "Ignore differences from Redux consensus and unmapping routines");

        configBuilder.addFlag(IGNORE_DUP_DIFFS, "Ignore duplicate diffs");
        configBuilder.addFlag(IGNORE_ALTERATIONS, "Ignore Redux consensus reads and unmappings");
        configBuilder.addFlag(IGNORE_CONSENSUS_READS, "Ignore consensus reads");
        configBuilder.addFlag(IGNORE_SUPPLEMENTARY_READS, "Ignore supplementary reads");
        configBuilder.addFlag(IGNORE_SECONDARY_READS, "Ignore secondary reads");
        configBuilder.addFlag(IGNORE_SUPPLEMENTARY_ATTRIBUTE, "Ignore supplementary attribute, can change from unmapping");
        configBuilder.addFlag(IGNORE_REDUX_UNMAPPED, "Ignore differences in reads unmapped by Redux");
        configBuilder.addFlag(CHECK_BASES_QUALS, "Ignore differences in consensus read bases and quals");
        configBuilder.addFlag(STANDARD_CHROMOSOMES, "Only process standard human chromosomes");

        configBuilder.addInteger(CIGAR_BOUNDARY_TOLERANCE,
                "Tolerate up to N bp drift on M-block boundaries when comparing CIGARs (0 = strict)", 0);

        configBuilder.addConfigItem(IGNORE_REGIONS_FILES,
                "Comma-separated TSV path(s) with Chromosome/PosStart/PosEnd columns; readnames whose any"
                        + " alignment overlaps these regions by more than half its reference span are dropped"
                        + " from comparison in both BAMs.");

        addRefGenomeFile(configBuilder, false);
        addSpecificChromosomesRegionsConfig(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
    }

    public boolean ignoreUnmapped() { return SpecificChrRegions != null && !SpecificChrRegions.Regions.isEmpty(); }

    @VisibleForTesting
    public CompareConfig(
            boolean ignoreDupDiffs, boolean ignoreConsensusReads, boolean ignoreSupplementaryReads,
            boolean ignoreSupplementaryAttribute, boolean ignoreReduxUnmapped, boolean ignoreReduxAlterations, boolean checkBaseAndQuals)
    {
        IgnoreDupDiffs = ignoreDupDiffs;
        IgnoreConsensusReads = ignoreConsensusReads;
        IgnoreSupplementaryReads = ignoreSupplementaryReads;
        IgnoreSupplementaryAttribute = ignoreSupplementaryAttribute;
        IgnoreReduxUnmapped = ignoreReduxUnmapped;
        IgnoreReduxAlterations = ignoreReduxAlterations;
        CheckBasesAndQuals = checkBaseAndQuals;
        IgnoreSecondaryReads = ignoreSupplementaryReads;

        CompareCoordsOnly = false;
        StandardChromosomes = false;
        CigarBoundaryTolerance = 0;

        OutputFile = null;
        OrigBamFile = null;
        NewBamFile = null;
        RefGenomeFile = null;
        RefGenVersion = null;
        PartitionSize = 0;
        MaxCachedReadsPerThread = 0;
        ExcludeRegions = false;
        Threads = 0;
        LogReadIds = null;
        SpecificChrRegions = null;
        IgnoreRegionFiles = Collections.emptyList();
    }
}
