package com.hartwig.hmftools.esvee.prep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.KNOWN_HOTSPOT_FILE;
import static com.hartwig.hmftools.esvee.common.FileCommon.PREP_FILE_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.formOutputFile;
import static com.hartwig.hmftools.esvee.common.FileCommon.parseBamFiles;
import static com.hartwig.hmftools.esvee.common.FileCommon.parseSampleIds;
import static com.hartwig.hmftools.esvee.common.FileCommon.registerCommonConfig;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SUPPORTING_READ_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_DISC_STATS_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_FRAG_LENGTH_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_JUNCTION_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.PREP_BAM;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.PREP_READ;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.parseConfigStr;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;
import com.hartwig.hmftools.esvee.prep.types.WriteType;

import htsjdk.samtools.ValidationStringency;

public class PrepConfig
{
    public final List<String> SampleIds;
    public final List<String> BamFiles;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final ReadFilters ReadFiltering;
    public final HotspotCache Hotspots;
    public final BlacklistLocations Blacklist;

    public final int PartitionSize;
    public final ValidationStringency BamStringency;

    public final String OutputDir;
    public final String OutputId;
    public final Set<WriteType> WriteTypes;

    public final String BamToolPath;

    public final int Threads;
    public final boolean UseCacheBam;
    public final boolean TrimReadId;

    // debug
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final boolean TrackRemotes;
    public final boolean PerfDebug;
    public final int MaxFragmentLengthOverride;

    public final boolean NoCleanUp;

    private boolean mIsValid;

    // both of these are set from the fragment length distribution
    private int mReadLength;
    private boolean mUnpairedReads;

    // config strings
    public static final String BAM_FILE = "bam_file";
    public static final String BAM_FILE_DESC = "BAM file paths separated by ','";
    public static final String SAMPLE_ID_DESC = "List of samples separated by ','";

    public static final String BLACKLIST_BED = "blacklist_bed";

    private static final String WRITE_TYPES = "write_types";

    public static final String READ_LENGTH = "read_length";
    public static final String PARTITION_SIZE = "partition_size";

    private static final String TRACK_REMOTES = "track_remotes";
    private static final String NO_CACHE_BAM = "no_cache_bam";
    private static final String NO_CLEAN_UP = "no_clean_up";
    private static final String NO_TRIM_READ_ID = "no_trim_read_id";
    private static final String MAX_FRAG_LENGTH_OVERRIDE = "max_frag_length_override";

    public PrepConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        SampleIds = parseSampleIds(configBuilder);
        BamFiles = parseBamFiles(configBuilder);

        if(SampleIds.isEmpty() || SampleIds.size() != BamFiles.size())
        {
            SV_LOGGER.error("samples({}) and bam({}) not matched or empty", SampleIds.size(), BamFiles.size());
            System.exit(1);
        }

        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
        }
        else
        {
            OutputDir = pathFromFile(BamFiles.get(0));
        }

        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(OutputDir == null || RefGenomeFile == null)
        {
            SV_LOGGER.error("missing config: refGenome({}) outputDir({})",
                    RefGenomeFile != null, OutputDir != null);
            System.exit(1);
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        SV_LOGGER.info("output({}) {}",
                OutputDir, OutputId != null ? format("outputId(%s)", OutputId) : "");

        Hotspots = new HotspotCache(configBuilder.getValue(KNOWN_HOTSPOT_FILE));
        Blacklist = new BlacklistLocations(configBuilder.getValue(BLACKLIST_BED));

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        ReadFiltering = new ReadFilters(ReadFilterConfig.from(configBuilder));

        WriteTypes = Sets.newHashSet(parseConfigStr(configBuilder.getValue(WRITE_TYPES)));

        BamStringency = BamUtils.validationStringency(configBuilder);
        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;

        LogReadIds = parseLogReadIds(configBuilder);

        Threads = parseThreads(configBuilder);

        MaxFragmentLengthOverride = configBuilder.getInteger(MAX_FRAG_LENGTH_OVERRIDE);

        // optimisations and debug
        TrimReadId = !configBuilder.hasFlag(NO_TRIM_READ_ID) && !SpecificChrRegions.hasFilters() && BamFiles.size() == 1;
        UseCacheBam = !configBuilder.hasFlag(NO_CACHE_BAM) && !SpecificChrRegions.hasFilters();
        TrackRemotes = configBuilder.hasFlag(TRACK_REMOTES);
        NoCleanUp = configBuilder.hasFlag(NO_CLEAN_UP);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);

        mReadLength = configBuilder.getInteger(READ_LENGTH);
        mUnpairedReads = false;
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        if(!Hotspots.isValid() || !Blacklist.isValid())
            return false;

        return true;
    }

    public int readLength() { return mReadLength; }
    public void setReadLength(int length) { mReadLength = length; }

    public boolean unpairedReads() { return mUnpairedReads; }
    public void setUnpairedReads(boolean unpaired) { mUnpairedReads = unpaired; }

    public String sampleId() { return SampleIds.get(0); }
    public String bamFile() { return BamFiles.get(0); }

    public String formFilename(final WriteType writeType) { return formFilename(writeType, sampleId()); }

    public String formFilename(final WriteType writeType, final String sampleId)
    {
        String fileExtension = "";

        switch(writeType)
        {
            case PREP_READ:
                fileExtension = "read" + TSV_EXTENSION;
                break;

            case UNSORTED_BAM:
                fileExtension = "unsorted.bam";
                break;

            case PREP_BAM:
                fileExtension = "bam";
                break;

            case CACHE_BAM:
                fileExtension = "cache";
                break;

            case PREP_JUNCTION:
                fileExtension = PREP_JUNCTION_FILE_ID;
                break;

            case FRAGMENT_LENGTH_DIST:
                fileExtension = PREP_FRAG_LENGTH_FILE_ID;
                break;

            case DISCORDANT_STATS:
                fileExtension = PREP_DISC_STATS_FILE_ID;
                break;
        }

        return formOutputFile(OutputDir, sampleId, PREP_FILE_ID, fileExtension, OutputId);
    }

    public boolean writeReads() { return WriteTypes.contains(PREP_BAM) || WriteTypes.contains(PREP_READ); }

    public PrepConfig(int partitionSize)
    {
        mIsValid = true;
        SampleIds = Collections.emptyList();
        BamFiles = Collections.emptyList();
        RefGenomeFile = "";
        OutputDir = null;
        OutputId = null;

        RefGenVersion = V37;

        Hotspots = new HotspotCache(null);
        Blacklist = new BlacklistLocations(null);

        PartitionSize = partitionSize;

        ReadFiltering = new ReadFilters(new ReadFilterConfig(
                MIN_ALIGNMENT_BASES,
                MIN_MAP_QUALITY,
                MIN_INSERT_ALIGNMENT_OVERLAP,
                MIN_SOFT_CLIP_LENGTH,
                MIN_SOFT_CLIP_HIGH_QUAL_PERC,
                MIN_SUPPORTING_READ_DISTANCE,
                MIN_INDEL_LENGTH,
                MIN_JUNCTION_SUPPORT));

        BamStringency = ValidationStringency.STRICT;
        WriteTypes = Sets.newHashSet();
        SpecificChrRegions = new SpecificRegions();
        LogReadIds = Lists.newArrayList();
        BamToolPath = null;
        Threads = 1;
        TrackRemotes = true;
        UseCacheBam = false;
        PerfDebug = false;
        TrimReadId = false;
        NoCleanUp = false;
        MaxFragmentLengthOverride = -1;

        mUnpairedReads = false;
        mReadLength = DEFAULT_READ_LENGTH;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_ID_DESC); // not required since may take from tumor and ref config values
        configBuilder.addPaths(BAM_FILE, false, BAM_FILE_DESC);

        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(KNOWN_HOTSPOT_FILE, false, "Known fusion hotspot BED file");
        configBuilder.addPath(BLACKLIST_BED, false, "Blacklist regions BED file");
        configBuilder.addInteger(READ_LENGTH, "Read length", DEFAULT_READ_LENGTH);
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);

        registerCommonConfig(configBuilder);

        configBuilder.addConfigItem(WRITE_TYPES, enumValueSelectionAsStr(WriteType.values(), "Write types"));

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        configBuilder.addFlag(NO_CACHE_BAM, "Write a BAM to cache candidate reads");
        configBuilder.addFlag(TRACK_REMOTES, "Track support for remote junctions");
        configBuilder.addFlag(NO_TRIM_READ_ID, "Disable use of a shortened readId internally");
        configBuilder.addFlag(NO_CLEAN_UP, "Keep candidate cache BAM files");
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addInteger(MAX_FRAG_LENGTH_OVERRIDE, "Set max fragment length instead of calculating", -1);
        addValidationStringencyOption(configBuilder);
        ReadFilterConfig.addConfig(configBuilder);
        BamToolName.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder, false);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
