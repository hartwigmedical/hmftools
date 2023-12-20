package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
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

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.samtools.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadFilterConfig;
import com.hartwig.hmftools.svprep.reads.ReadFilters;

import org.apache.commons.cli.ParseException;

import htsjdk.samtools.ValidationStringency;

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
    public int ReadLength; // can be set from default, config or the fragment length distribution routine
    public final boolean CalcFragmentLength;
    public final ValidationStringency BamStringency;

    public final String OutputDir;
    public final String OutputId;
    public final Set<WriteType> WriteTypes;

    public final int Threads;
    public final boolean UseCacheBam;
    public final boolean TrimReadId;
    public final boolean UnpairedReads;

    // debug
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final boolean TrackRemotes;
    public final boolean PerfDebug;

    // throttling and down-sampling - off by default
    public final int JunctionFragmentCap;
    public final int MaxPartitionReads;
    public final boolean CaptureDepth;
    public final boolean NoCleanUp;

    private boolean mIsValid;

    // config strings
    public static final String BAM_FILE = "bam_file";
    private static final String KNOWN_FUSION_BED = "known_fusion_bed";
    public static final String BLACKLIST_BED = "blacklist_bed";
    private static final String EXISTING_JUNCTION_FILE = "existing_junction_file";

    private static final String WRITE_TYPES = "write_types";

    public static final String READ_LENGTH = "read_length";
    private static final String CALC_FRAG_LENGTH = "calc_fragment_length";
    private static final String PARTITION_SIZE = "partition_size";

    private static final String MAX_PARTITION_READS = "max_partition_reads";
    private static final String CAPTURE_DEPTH = "capture_depth";
    private static final String TRACK_REMOTES = "track_remotes";
    private static final String NO_CACHE_BAM = "no_cache_bam";
    private static final String NO_CLEAN_UP = "no_clean_up";
    private static final String NO_TRIM_READ_ID = "no_trim_read_id";
    private static final String JUNCTION_FRAGS_CAP = "junction_frags_cap";
    private static final String UNPAIRED_READS = "unpaired_reads";

    public SvConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        SampleId = configBuilder.getValue(SAMPLE);
        BamFile = configBuilder.getValue(BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
        }
        else
        {
            OutputDir = pathFromFile(BamFile);
        }

        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(SampleId == null || BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            SV_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    SampleId != null, BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        SV_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        SV_LOGGER.info("output({})", OutputDir);

        Hotspots = new HotspotCache(configBuilder.getValue(KNOWN_FUSION_BED));
        Blacklist = new BlacklistLocations(configBuilder.getValue(BLACKLIST_BED));

        ExistingJunctionFile = configBuilder.getValue(EXISTING_JUNCTION_FILE);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        ReadLength = configBuilder.getInteger(READ_LENGTH);

        ReadFiltering = new ReadFilters(ReadFilterConfig.from(configBuilder));

        WriteTypes = Sets.newHashSet();

        if(configBuilder.hasValue(WRITE_TYPES))
        {
            String[] writeTypes = configBuilder.getValue(WRITE_TYPES).split(ITEM_DELIM, -1);
            Arrays.stream(writeTypes).forEach(x -> WriteTypes.add(WriteType.valueOf(x)));
        }
        else
        {
            WriteTypes.add(WriteType.JUNCTIONS);
            WriteTypes.add(WriteType.BAM);
        }

        CalcFragmentLength = configBuilder.hasFlag(CALC_FRAG_LENGTH) || WriteTypes.contains(FRAGMENT_LENGTH_DIST);
        BamStringency = BamUtils.validationStringency(configBuilder);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;

        LogReadIds = parseLogReadIds(configBuilder);

        Threads = parseThreads(configBuilder);

        // optimisations and debug
        TrimReadId = !configBuilder.hasFlag(NO_TRIM_READ_ID) && !SpecificChrRegions.hasFilters();
        UnpairedReads = configBuilder.hasFlag(UNPAIRED_READS);
        UseCacheBam = !configBuilder.hasFlag(NO_CACHE_BAM) && !SpecificChrRegions.hasFilters();
        MaxPartitionReads = configBuilder.getInteger(MAX_PARTITION_READS);
        JunctionFragmentCap = configBuilder.getInteger(JUNCTION_FRAGS_CAP);
        CaptureDepth = configBuilder.hasFlag(CAPTURE_DEPTH);
        TrackRemotes = configBuilder.hasFlag(TRACK_REMOTES);
        NoCleanUp = configBuilder.hasFlag(NO_CLEAN_UP);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        if(!Hotspots.isValid() || !Blacklist.isValid())
            return false;

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
            case READS: return filename + "reads" + TSV_EXTENSION;
            case BAM: return filename + "bam";
            case CACHE_BAM: return filename + "cache";
            case JUNCTIONS: return filename + "junctions" + TSV_EXTENSION;
            case FRAGMENT_LENGTH_DIST: return filename + "fragment_lengths" + TSV_EXTENSION;
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
        BamStringency = ValidationStringency.STRICT;
        WriteTypes = Sets.newHashSet();
        SpecificChrRegions = new SpecificRegions();
        LogReadIds = Lists.newArrayList();
        Threads = 1;
        MaxPartitionReads = 0;
        TrackRemotes = true;
        UseCacheBam = false;
        PerfDebug = false;
        TrimReadId = false;
        UnpairedReads = false;
        NoCleanUp = false;
        JunctionFragmentCap = 0;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, true, "BAM file location");
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(KNOWN_FUSION_BED, false, "Known fusion hotspot BED file");
        configBuilder.addPath(BLACKLIST_BED, false, "Blacklist regions BED file");
        configBuilder.addPath(EXISTING_JUNCTION_FILE, false, "Load existing junction file to find supporting reads");
        configBuilder.addInteger(READ_LENGTH, "Read length", DEFAULT_READ_LENGTH);
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addFlag(CALC_FRAG_LENGTH, "Calculate distribution for fragment length");
        configBuilder.addConfigItem(WRITE_TYPES, "Write types: " + WriteType.values().toString());
        configBuilder.addFlag(UNPAIRED_READS, "Unpaired reads ignores non-expect junction support");
        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        configBuilder.addInteger(MAX_PARTITION_READS, "Limit to stop processing reads in partition, for debug", 0);
        configBuilder.addFlag(CAPTURE_DEPTH, "Capture depth for junctions");
        configBuilder.addFlag(NO_CACHE_BAM, "Write a BAM to cache candidate reads");
        configBuilder.addFlag(TRACK_REMOTES, "Track support for remote junctions");
        configBuilder.addFlag(NO_TRIM_READ_ID, "Disable use of a shortened readId internally");
        configBuilder.addFlag(NO_CLEAN_UP, "Keep candidate cache BAM files");
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addInteger(JUNCTION_FRAGS_CAP, "Limit to supporting reads added to a junction", 0);
        addValidationStringencyOption(configBuilder);
        ReadFilterConfig.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder, false);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
