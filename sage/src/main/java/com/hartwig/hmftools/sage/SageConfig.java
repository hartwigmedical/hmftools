package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.sage.SageCommon.SAMPLE_DELIM;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_PARTITION_SLICES;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH_PANEL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_SLICE_SIZE;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.samtools.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.quality.QualityConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationConfig;

import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.ValidationStringency;

public class SageConfig
{
    public final List<String> ReferenceIds;
    public final List<String> ReferenceBams;

    public final String SampleDataDir;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String OutputFile;
    public final FilterConfig Filter;
    public final QualityConfig Quality;
    public final QualityRecalibrationConfig QualityRecalibration;
    public final boolean IncludeMT;
    public final boolean SyncFragments;
    public final int RegionSliceSize;
    public final int MinMapQuality;
    public final int MaxReadDepth;
    public final int MaxReadDepthPanel;
    public final int ReadContextFlankSize;
    public final int MaxPartitionSlices;
    public final ValidationStringency BamStringency;

    public final String Version;
    public final int Threads;

    public final boolean TrackUMIs;

    // debug
    public final SpecificRegions SpecificChrRegions;
    public final Set<Integer> SpecificPositions;
    public final boolean LogEvidenceReads;
    public final boolean LogLpsData;
    public final double PerfWarnTime;

    private boolean mIsValid;

    // allow Sage to auto-adjust to the longest observed read length
    private int mReadLength;

    private static final String REFERENCE = "reference";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String OUTPUT_VCF = "out";
    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String MAX_READ_DEPTH = "max_read_depth";
    private static final String MAX_READ_DEPTH_PANEL = "max_read_depth_panel";
    private static final String SLICE_SIZE = "slice_size";
    private static final String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";
    private static final String INCLUDE_MT = "include_mt";
    private static final String READ_LENGTH = "read_length";
    private static final String SYNC_FRAGMENTS = "sync_fragments";
    private static final String TRACK_UMIS = "track_umis";
    private static final String MAX_PARTITION_SLICES = "max_partition_slices";

    private static final String SPECIFIC_POSITIONS = "specific_positions";
    private static final String LOG_EVIDENCE_READS = "log_evidence_reads";
    private static final String LOG_LPS_DATA = "log_lps_data";
    private static final String PERF_WARN_TIME = "perf_warn_time";

    public SageConfig(final String version, final ConfigBuilder configBuilder)
    {
        mIsValid = true;
        Version = version;
        
        RefGenVersion = RefGenomeVersion.from(configBuilder);

        ReferenceIds = Lists.newArrayList();
        if(configBuilder.hasValue(REFERENCE))
        {
            ReferenceIds.addAll(Arrays.asList(configBuilder.getValue(REFERENCE).split(SAMPLE_DELIM)));
        }

        SampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG, ""));

        ReferenceBams = Lists.newArrayList();

        if(configBuilder.hasValue(REFERENCE_BAM))
        {
            Arrays.stream(configBuilder.getValue(REFERENCE_BAM, Strings.EMPTY).split(SAMPLE_DELIM))
                    .forEach(x -> ReferenceBams.add(SampleDataDir + x));
        }


        IncludeMT = configBuilder.hasFlag(INCLUDE_MT);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;

        OutputFile = SampleDataDir + configBuilder.getValue(OUTPUT_VCF);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        BamStringency = BamUtils.validationStringency(configBuilder);
        RegionSliceSize = configBuilder.getInteger(SLICE_SIZE);
        ReadContextFlankSize = configBuilder.getInteger(READ_CONTEXT_FLANK_SIZE);
        MinMapQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxReadDepth = configBuilder.getInteger(MAX_READ_DEPTH);
        MaxReadDepthPanel = configBuilder.getInteger(MAX_READ_DEPTH_PANEL);

        mReadLength = configBuilder.getInteger(READ_LENGTH);

        MaxPartitionSlices = configBuilder.getInteger(MAX_PARTITION_SLICES);
        SyncFragments = configBuilder.hasFlag(SYNC_FRAGMENTS);

        Filter = new FilterConfig(configBuilder);
        Quality = new QualityConfig(configBuilder);
        QualityRecalibration = new QualityRecalibrationConfig(configBuilder);

        TrackUMIs = configBuilder.hasFlag(TRACK_UMIS);

        SpecificPositions = Sets.newHashSet();
        if(configBuilder.hasValue(SPECIFIC_POSITIONS))
        {
            final String positionList = configBuilder.getValue(SPECIFIC_POSITIONS, Strings.EMPTY);
            if(!positionList.isEmpty())
            {
                Arrays.stream(positionList.split(ITEM_DELIM)).forEach(x -> SpecificPositions.add(Integer.parseInt(x)));
            }
        }

        LogLpsData = configBuilder.hasFlag(LOG_LPS_DATA);
        LogEvidenceReads = SpecificChrRegions.hasFilters() && configBuilder.hasFlag(LOG_EVIDENCE_READS);

        if(LogEvidenceReads)
        {
            SG_LOGGER.trace("READ_EV,SampleId,Chromosome,Position,Ref,Alt,MatchType,ReadId,ReadStart,Cigar,LeftCore,Index,RightCore,ReadIndex,Quality");
        }

        PerfWarnTime = configBuilder.getDecimal(PERF_WARN_TIME);

        Threads = parseThreads(configBuilder);
    }

    public int getReadLength() { return mReadLength; }

    public void setReadLength(int readLength)
    {
        if(readLength != DEFAULT_READ_LENGTH)
        {
            SG_LOGGER.info("max observed read length set({})", readLength);
        }

        mReadLength = readLength;
    }

    public String formOutputDir()
    {
        return OutputFile != null ? new File(OutputFile).getParent() + File.separator : "";
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        if(ReferenceIds.size() != ReferenceBams.size())
        {
            SG_LOGGER.error("Each reference sample must have matching bam");
            return false;
        }

        for(String referenceBam : ReferenceBams)
        {
            if(!new File(referenceBam).exists())
            {
                SG_LOGGER.error("Unable to locate reference bam({})", referenceBam);
                return false;
            }
        }

        if(RefGenomeFile.isEmpty())
        {
            SG_LOGGER.error("Reference genome required");
            return false;
        }

        if(OutputFile.isEmpty())
        {
            SG_LOGGER.error("No output VCF file specified");
            return false;
        }

        final File outputDir = new File(OutputFile).getParentFile();
        if(outputDir != null && !outputDir.exists() && !outputDir.mkdirs())
        {
            SG_LOGGER.error("unable to write directory({})", outputDir.toString());
            return false;
        }

        return true;
    }

    public boolean processChromosome(final String chromosome)
    {
        if(SpecificChrRegions.excludeChromosome(chromosome))
            return false;

        if(HumanChromosome.contains(chromosome))
            return true;

        if(IncludeMT && MitochondrialChromosome.contains(chromosome))
            return true;

        return false;
    }

    public boolean logPerfStats() { return PerfWarnTime > 0; }

    public static void registerCommonConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REFERENCE, false, "Reference sample, or collection separated by ','");
        configBuilder.addConfigItem(REFERENCE_BAM, false, "Reference bam file");

        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, "Path to sample data files");
        configBuilder.addConfigItem(OUTPUT_VCF, true, "Output vcf");

        addRefGenomeConfig(configBuilder, true);

        // is this common?
        configBuilder.addInteger(
                READ_CONTEXT_FLANK_SIZE, "Size of read context flank", DEFAULT_READ_CONTEXT_FLANK_SIZE);

        configBuilder.addInteger(MIN_MAP_QUALITY, "Min map quality to apply to non-hotspot variants", DEFAULT_MIN_MAP_QUALITY);
        configBuilder.addInteger(READ_LENGTH, "Read length, otherwise will sample from BAM", 0);
        configBuilder.addFlag(INCLUDE_MT, "Call MT variants");
        configBuilder.addInteger(SLICE_SIZE, "Slice size", DEFAULT_SLICE_SIZE);
        configBuilder.addInteger(MAX_PARTITION_SLICES, "Max slices per partition", DEFAULT_MAX_PARTITION_SLICES);

        configBuilder.addInteger(MAX_READ_DEPTH, "Max depth to look for evidence", DEFAULT_MAX_READ_DEPTH);
        configBuilder.addInteger(MAX_READ_DEPTH_PANEL, "Max depth to look for evidence in panel", DEFAULT_MAX_READ_DEPTH_PANEL);
        configBuilder.addFlag(SYNC_FRAGMENTS, "Handle overlapping fragment reads in evidence phase");
        configBuilder.addFlag(TRACK_UMIS, "Record counts of UMI types");
        addValidationStringencyOption(configBuilder);

        FilterConfig.registerConfig(configBuilder);
        QualityConfig.registerConfig(configBuilder);
        QualityRecalibrationConfig.registerConfig(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addFlag(LOG_EVIDENCE_READS, "Write evidence read data");

        configBuilder.addFlag(LOG_LPS_DATA, "Log local phasing data");
        configBuilder.addDecimal(PERF_WARN_TIME, "Log details of partitions taking longer than X seconds", 0.0);

        // debug
        configBuilder.addConfigItem(SPECIFIC_POSITIONS, "Run for specific positions(s) separated by ';', for debug purposes");

        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    public SageConfig()
    {
        SampleDataDir = "";
        ReferenceIds = Lists.newArrayList();
        ReferenceBams = Lists.newArrayList();
        Filter = new FilterConfig();
        Quality = new QualityConfig();
        QualityRecalibration = new QualityRecalibrationConfig();
        SpecificChrRegions = new SpecificRegions();
        IncludeMT = false;
        RegionSliceSize = 500_000;
        MinMapQuality = DEFAULT_MIN_MAP_QUALITY;
        MaxReadDepth = DEFAULT_MAX_READ_DEPTH;
        MaxReadDepthPanel = DEFAULT_MAX_READ_DEPTH_PANEL;
        ReadContextFlankSize = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        mReadLength = DEFAULT_READ_LENGTH;
        MaxPartitionSlices = 1;
        RefGenomeFile = "refGenome";
        OutputFile = "out.vcf";
        Version = "1.0";
        Threads = 1;
        LogLpsData = false;
        PerfWarnTime = 0;
        RefGenVersion = V37;
        BamStringency = ValidationStringency.DEFAULT_STRINGENCY;
        TrackUMIs = false;
        SyncFragments = false;
        SpecificPositions = Sets.newHashSet();
        LogEvidenceReads = false;
    }
}
