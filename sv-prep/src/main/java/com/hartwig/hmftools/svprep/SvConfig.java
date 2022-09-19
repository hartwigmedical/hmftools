package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.svprep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SUPPORTING_READ_DISTANCE;
import static com.hartwig.hmftools.svprep.WriteType.BAM;
import static com.hartwig.hmftools.svprep.WriteType.FRAGMENT_LENGTH_DIST;
import static com.hartwig.hmftools.svprep.WriteType.READS;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadFilterConfig;
import com.hartwig.hmftools.svprep.reads.ReadFilters;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class SvConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final ReadFilters ReadFiltering;
    public final HotspotCache Hotspots;
    public final BlacklistLocations Blacklist;
    public final String ExistingJunctionFile;

    public final int PartitionSize;
    public final int ReadLength;
    public final boolean CalcFragmentLength;

    public final String OutputDir;
    public final String OutputId;
    public final Set<WriteType> WriteTypes;

    public final int Threads;
    public final boolean UseCacheBam;
    public final boolean TrimReadId;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean TrackRemotes;
    public final boolean PerfDebug;

    // throttling and down-sampling - off by default
    public final int JunctionFragmentCap;
    public final int MaxPartitionReads;
    public final boolean ApplyDownsampling;
    public final boolean CaptureDepth;
    public final boolean NoCleanUp;

    private boolean mIsValid;

    // config strings
    public static final String SAMPLE = "sample";
    private static final String BAM_FILE = "bam_file";
    private static final String KNOWN_FUSION_BED = "known_fusion_bed";
    public static final String BLACKLIST_BED = "blacklist_bed";
    private static final String EXISTING_JUNCTION_FILE = "existing_junction_file";

    private static final String WRITE_TYPES = "write_types";

    private static final String READ_LENGTH = "read_length";
    private static final String CALC_FRAG_LENGTH = "calc_fragment_length";
    private static final String PARTITION_SIZE = "partition_size";

    private static final String LOG_READ_IDS = "log_read_ids";
    private static final String MAX_PARTITION_READS = "max_partition_reads";
    private static final String APPLY_DOWNSAMPLING = "apply_downsampling";
    private static final String CAPTURE_DEPTH = "capture_depth";
    private static final String TRACK_REMOTES = "track_remotes";
    private static final String NO_CACHE_BAM = "no_cache_bam";
    private static final String NO_CLEAN_UP = "no_clean_up";
    private static final String NO_TRIM_READ_ID = "no_trim_read_id";
    private static final String PERF_DEBUG = "perf_debug";
    private static final String JUNCTION_FRAGS_CAP = "junction_frags_cap";

    public SvConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        if(SampleId == null || BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            SV_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    SampleId != null, BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        SV_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        SV_LOGGER.info("output({})", OutputDir);

        Hotspots = new HotspotCache(cmd.getOptionValue(KNOWN_FUSION_BED));
        Blacklist = new BlacklistLocations(cmd.getOptionValue(BLACKLIST_BED));

        ExistingJunctionFile = cmd.getOptionValue(EXISTING_JUNCTION_FILE);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));
        ReadLength = Integer.parseInt(cmd.getOptionValue(READ_LENGTH, String.valueOf(DEFAULT_READ_LENGTH)));

        ReadFiltering = new ReadFilters(ReadFilterConfig.from(cmd));

        WriteTypes = Sets.newHashSet();

        if(cmd.hasOption(WRITE_TYPES))
        {
            String[] writeTypes = cmd.getOptionValue(WRITE_TYPES).split(ITEM_DELIM, -1);
            Arrays.stream(writeTypes).forEach(x -> WriteTypes.add(WriteType.valueOf(x)));
        }
        else
        {
            WriteTypes.add(WriteType.JUNCTIONS);
            WriteTypes.add(WriteType.BAM);
        }

        CalcFragmentLength = cmd.hasOption(CALC_FRAG_LENGTH) || WriteTypes.contains(FRAGMENT_LENGTH_DIST);

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        try
        {
            loadSpecificChromsomesOrRegions(cmd, SpecificChromosomes, SpecificRegions, SV_LOGGER);
        }
        catch(ParseException e)
        {
            mIsValid = false;
        }

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();

        Threads = parseThreads(cmd);

        // optimisations and debug
        TrimReadId = !cmd.hasOption(NO_TRIM_READ_ID) && SpecificRegions.isEmpty();
        UseCacheBam = !cmd.hasOption(NO_CACHE_BAM) && SpecificRegions.isEmpty();
        MaxPartitionReads = Integer.parseInt(cmd.getOptionValue(MAX_PARTITION_READS, "0"));
        JunctionFragmentCap = Integer.parseInt(cmd.getOptionValue(JUNCTION_FRAGS_CAP, "0"));
        CaptureDepth = cmd.hasOption(CAPTURE_DEPTH);
        ApplyDownsampling = cmd.hasOption(APPLY_DOWNSAMPLING);
        TrackRemotes = cmd.hasOption(TRACK_REMOTES);
        NoCleanUp = cmd.hasOption(NO_CLEAN_UP);
        PerfDebug = cmd.hasOption(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        if(!Files.exists(Paths.get(BamFile)))
        {
            SV_LOGGER.error("invalid bam file path: {}", BamFile);
            return false;
        }

        if(!Files.exists(Paths.get(RefGenomeFile)))
        {
            SV_LOGGER.error("invalid ref genome file: {}", RefGenomeFile);
            return false;
        }

        if(!Hotspots.isValid() || !Blacklist.isValid())
            return false;

        if(ExistingJunctionFile != null && !Files.exists(Paths.get(ExistingJunctionFile)))
        {
            SV_LOGGER.error("invalid existing junctions file: {}", ExistingJunctionFile);
            return false;
        }

        return true;
    }

    public String formFilename(final WriteType writeType)
    {
        String filename = OutputDir + SampleId;

        filename += ".sv_prep.";

        if(OutputId != null)
            filename += OutputId + ".";

        switch(writeType)
        {
            case READS: return filename + "reads.csv";
            case BAM: return filename + "bam";
            case CACHE_BAM: return filename + "cache";
            case JUNCTIONS: return filename + "junctions.csv";
            case FRAGMENT_LENGTH_DIST: return filename + "fragment_lengths";
        }

        return null;
    }

    public boolean writeReads() { return WriteTypes.contains(BAM) || WriteTypes.contains(READS); }

    public SvConfig(int partitionSize)
    {
        mIsValid = true;
        SampleId = "TEST";
        BamFile = null;
        RefGenomeFile = "";
        OutputDir = null;
        OutputId = null;

        RefGenVersion = V37;

        Hotspots = new HotspotCache(null);
        Blacklist = new BlacklistLocations(null);
        ExistingJunctionFile = null;

        PartitionSize = partitionSize;

        ReadLength = DEFAULT_READ_LENGTH;

        ReadFiltering = new ReadFilters(new ReadFilterConfig(
                MIN_ALIGNMENT_BASES,
                MIN_MAP_QUALITY,
                MIN_INSERT_ALIGNMENT_OVERLAP,
                MIN_SOFT_CLIP_LENGTH,
                MIN_SOFT_CLIP_MIN_BASE_QUAL,
                MIN_SOFT_CLIP_HIGH_QUAL_PERC,
                MIN_SUPPORTING_READ_DISTANCE,
                MIN_INDEL_LENGTH,
                MIN_JUNCTION_SUPPORT));

        CalcFragmentLength = false;
        CaptureDepth = false;
        WriteTypes = Sets.newHashSet();
        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();
        LogReadIds = Lists.newArrayList();
        Threads = 1;
        MaxPartitionReads = 0;
        ApplyDownsampling = false;
        TrackRemotes = true;
        UseCacheBam = false;
        PerfDebug = false;
        TrimReadId = false;
        NoCleanUp = false;
        JunctionFragmentCap = 0;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "RNA BAM file location");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default) or 38");
        options.addOption(KNOWN_FUSION_BED, true, "Known fusion hotspot BED file");
        options.addOption(BLACKLIST_BED, true, "Blacklist regions BED file");
        options.addOption(EXISTING_JUNCTION_FILE, true, "Load existing junction file to find supporting reads");
        options.addOption(READ_LENGTH, true, "Read length, default: " + DEFAULT_READ_LENGTH);
        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_CHR_PARTITION_SIZE);
        options.addOption(CALC_FRAG_LENGTH, false, "Calculate distribution for fragment length");
        options.addOption(WRITE_TYPES, true, "Write types: " + WriteType.values().toString());
        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(MAX_PARTITION_READS, true, "Limit to stop processing reads in partition, for debug");
        options.addOption(CAPTURE_DEPTH, false, "Capture depth for junctions");
        options.addOption(APPLY_DOWNSAMPLING, false, "Apply downsampling of reads in high-depth regions");
        options.addOption(NO_CACHE_BAM, false, "Write a BAM to cache candidate reads");
        options.addOption(TRACK_REMOTES, false, "Track support for remote junctions");
        options.addOption(NO_TRIM_READ_ID, false, "Use a shortened readId internally");
        options.addOption(NO_CLEAN_UP, false, "Keep candidate cache files");
        options.addOption(PERF_DEBUG, false, "Detailed performance tracking and logging");
        options.addOption(JUNCTION_FRAGS_CAP, true, "Limit to supporting reads added to a junction");
        ReadFilterConfig.addCmdLineArgs(options);
        addThreadOptions(options);

        return options;
    }
}
