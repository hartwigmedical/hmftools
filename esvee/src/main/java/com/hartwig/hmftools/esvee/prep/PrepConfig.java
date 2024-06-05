package com.hartwig.hmftools.esvee.prep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.BamToolName.BAMTOOL_PATH;
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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.PREP_FILE_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.formOutputFile;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SUPPORTING_READ_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_FRAG_LENGTH_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_JUNCTIONS_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.BAM;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.FRAGMENT_LENGTH_DIST;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.READS;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;
import com.hartwig.hmftools.esvee.prep.types.ReadFilters;
import com.hartwig.hmftools.esvee.prep.types.WriteType;

import org.jetbrains.annotations.Nullable;

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
    public int ReadLength; // can be set from default, config or the fragment length distribution routine
    public final boolean CalcFragmentLength;
    public final ValidationStringency BamStringency;

    public final String OutputDir;
    public final String OutputId;
    public final Set<WriteType> WriteTypes;

    public final String BamToolPath;

    public final int Threads;
    public final boolean UseCacheBam;
    public final boolean TrimReadId;
    public final boolean UnpairedReads;

    // debug
    public final SpecificRegions SpecificChrRegions;
    public final List<String> LogReadIds;
    public final boolean TrackRemotes;
    public final boolean PerfDebug;

    public final boolean NoCleanUp;

    private boolean mIsValid;

    // config strings
    public static final String BAM_FILE = "bam_file";
    public static final String BAM_FILES = "bam_files";
    private static final String KNOWN_FUSION_BED = "known_fusion_bed";
    public static final String BLACKLIST_BED = "blacklist_bed";

    private static final String WRITE_TYPES = "write_types";

    public static final String READ_LENGTH = "read_length";
    private static final String CALC_FRAG_LENGTH = "calc_fragment_length";
    private static final String PARTITION_SIZE = "partition_size";

    private static final String TRACK_REMOTES = "track_remotes";
    private static final String NO_CACHE_BAM = "no_cache_bam";
    private static final String NO_CLEAN_UP = "no_clean_up";
    private static final String NO_TRIM_READ_ID = "no_trim_read_id";
    private static final String UNPAIRED_READS = "unpaired_reads";

    public PrepConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        SampleIds = Arrays.stream(configBuilder.getValue(SAMPLE).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());

        if(configBuilder.hasValue(BAM_FILE))
            BamFiles = Lists.newArrayList(configBuilder.getValue(BAM_FILE));
        else
            BamFiles = Arrays.stream(configBuilder.getValue(BAM_FILES).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());

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

        if(SampleIds.size() != BamFiles.size())
        {
            SV_LOGGER.error("samples({}) and bam({}) not matched", SampleIds.size(), BamFiles.size());
            mIsValid = false;
        }

        if(SampleIds.isEmpty() || BamFiles.isEmpty() || OutputDir == null || RefGenomeFile == null)
        {
            SV_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    !SampleIds.isEmpty(), !BamFiles.isEmpty(), RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        SV_LOGGER.info("output({}) {}",
                OutputDir, OutputId != null ? format("outputId(%s)", OutputId) : "");

        Hotspots = new HotspotCache(configBuilder.getValue(KNOWN_FUSION_BED));
        Blacklist = new BlacklistLocations(configBuilder.getValue(BLACKLIST_BED));

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
        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;

        LogReadIds = parseLogReadIds(configBuilder);

        Threads = parseThreads(configBuilder);

        // optimisations and debug
        TrimReadId = !configBuilder.hasFlag(NO_TRIM_READ_ID) && !SpecificChrRegions.hasFilters();
        UnpairedReads = configBuilder.hasFlag(UNPAIRED_READS);
        UseCacheBam = !configBuilder.hasFlag(NO_CACHE_BAM) && !SpecificChrRegions.hasFilters();
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

    public String sampleId() { return SampleIds.get(0); }
    public String bamFile() { return BamFiles.get(0); }

    public String formFilename(final WriteType writeType) { return formFilename(writeType, sampleId()); }

    public String formFilename(final WriteType writeType, final String sampleId)
    {
        String fileExtension = "";

        switch(writeType)
        {
            case READS:
                fileExtension = "reads" + TSV_EXTENSION;
                break;

            case UNSORTED_BAM:
                fileExtension = "unsorted.bam";
                break;

            case BAM:
                fileExtension = "bam";
                break;

            case CACHE_BAM:
                fileExtension = "cache";
                break;

            case JUNCTIONS:
                fileExtension = PREP_JUNCTIONS_FILE_ID;
                break;

            case FRAGMENT_LENGTH_DIST:
                fileExtension = PREP_FRAG_LENGTH_FILE_ID;
                break;
        }

        return formOutputFile(OutputDir, sampleId, PREP_FILE_ID, fileExtension, OutputId);
    }

    public boolean writeReads() { return WriteTypes.contains(BAM) || WriteTypes.contains(READS); }

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

        ReadLength = DEFAULT_READ_LENGTH;

        ReadFiltering = new ReadFilters(new ReadFilterConfig(
                MIN_ALIGNMENT_BASES,
                MIN_MAP_QUALITY,
                MIN_INSERT_ALIGNMENT_OVERLAP,
                MIN_SOFT_CLIP_LENGTH,
                LOW_BASE_QUAL_THRESHOLD,
                MIN_SOFT_CLIP_HIGH_QUAL_PERC,
                MIN_SUPPORTING_READ_DISTANCE,
                MIN_INDEL_LENGTH,
                MIN_JUNCTION_SUPPORT));

        CalcFragmentLength = false;
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
        UnpairedReads = false;
        NoCleanUp = false;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, false, "BAM file location");
        configBuilder.addPaths(BAM_FILES, false, "BAM file(s) location");
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(KNOWN_FUSION_BED, false, "Known fusion hotspot BED file");
        configBuilder.addPath(BLACKLIST_BED, false, "Blacklist regions BED file");
        configBuilder.addInteger(READ_LENGTH, "Read length", DEFAULT_READ_LENGTH);
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addFlag(CALC_FRAG_LENGTH, "Calculate distribution for fragment length");
        configBuilder.addConfigItem(WRITE_TYPES, "Write types: " + WriteType.values().toString());
        configBuilder.addFlag(UNPAIRED_READS, "Unpaired reads ignores non-expect junction support");
        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        configBuilder.addFlag(NO_CACHE_BAM, "Write a BAM to cache candidate reads");
        configBuilder.addFlag(TRACK_REMOTES, "Track support for remote junctions");
        configBuilder.addFlag(NO_TRIM_READ_ID, "Disable use of a shortened readId internally");
        configBuilder.addFlag(NO_CLEAN_UP, "Keep candidate cache BAM files");
        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        addValidationStringencyOption(configBuilder);
        ReadFilterConfig.addConfig(configBuilder);
        BamToolName.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder, false);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
