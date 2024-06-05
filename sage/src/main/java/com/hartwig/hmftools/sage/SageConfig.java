package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAMS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.sage.SageCommon.SAMPLE_DELIM;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_PARTITION_SLICES;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH_PANEL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_SLICE_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.VIS_VARIANT_BUFFER;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sage.bqr.BqrConfig;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.quality.QualityConfig;
import com.hartwig.hmftools.sage.vis.VisConfig;

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
    public final BqrConfig BQR;
    public final String JitterParamsDir;
    public final boolean IncludeMT;
    public final boolean SyncFragments;
    public final int RegionSliceSize;
    public final int MinMapQuality;
    public final int MaxReadDepth;
    public final int MaxReadDepthPanel;
    public final int ReadContextFlankLength;
    public final int MaxPartitionSlices;
    public final ValidationStringency BamStringency;
    public final SequencingConfig Sequencing;

    public final VisConfig Visualiser;

    public final String Version;
    public final int Threads;

    public final boolean WriteFragmentLengths;

    // debug
    public final SpecificRegions SpecificChrRegions;
    public final List<BasePosition> SpecificPositions;
    public final boolean LogEvidenceReads;
    public final boolean LogLpsData;
    public final double PerfWarnTime;

    private boolean mIsValid;

    // allow Sage to auto-adjust to the longest observed read length
    private int mReadLength;

    private static final String OUTPUT_VCF = "output_vcf";
    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String MAX_READ_DEPTH = "max_read_depth";
    private static final String MAX_READ_DEPTH_PANEL = "max_read_depth_panel";
    private static final String SLICE_SIZE = "slice_size";
    private static final String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";
    private static final String INCLUDE_MT = "include_mt";
    private static final String READ_LENGTH = "read_length";
    private static final String NO_FRAGMENT_SYNC = "no_fragment_sync";
    private static final String WRITE_FRAG_LENGTHS = "write_frag_lengths";
    private static final String MAX_PARTITION_SLICES = "max_partition_slices";
    private static final String JITTER_PARAMS_DIR = "jitter_param_dir";

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

        OutputFile = SampleDataDir + configBuilder.getValue(OUTPUT_VCF);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        BamStringency = BamUtils.validationStringency(configBuilder);
        RegionSliceSize = configBuilder.getInteger(SLICE_SIZE);
        ReadContextFlankLength = configBuilder.getInteger(READ_CONTEXT_FLANK_SIZE);

        MaxReadDepth = configBuilder.getInteger(MAX_READ_DEPTH);
        MaxReadDepthPanel = configBuilder.getInteger(MAX_READ_DEPTH_PANEL);

        mReadLength = configBuilder.getInteger(READ_LENGTH);

        MaxPartitionSlices = configBuilder.getInteger(MAX_PARTITION_SLICES);
        SyncFragments = !configBuilder.hasFlag(NO_FRAGMENT_SYNC);

        Filter = new FilterConfig(configBuilder);
        Quality = new QualityConfig(configBuilder);
        BQR = new BqrConfig(configBuilder);
        JitterParamsDir = configBuilder.getValue(JITTER_PARAMS_DIR);

        MinMapQuality = configBuilder.getInteger(MIN_MAP_QUALITY);

        Sequencing = SequencingConfig.from(configBuilder);

        WriteFragmentLengths = configBuilder.hasFlag(WRITE_FRAG_LENGTHS);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;

        Visualiser = new VisConfig(configBuilder, outputDir());

        if(Visualiser.Enabled && !Visualiser.SpecificVariants.isEmpty() && SpecificChrRegions.Regions.isEmpty())
        {
            for(SimpleVariant visVariant : Visualiser.SpecificVariants)
            {
                SpecificChrRegions.addRegion(new ChrBaseRegion(
                        visVariant.Chromosome, visVariant.Position - VIS_VARIANT_BUFFER,
                        visVariant.Position + VIS_VARIANT_BUFFER));
            }
        }

        SpecificPositions = Lists.newArrayList();

        if(configBuilder.hasValue(SPECIFIC_POSITIONS))
        {
            String[] specPositionsStr = configBuilder.getValue(SPECIFIC_POSITIONS).split(ITEM_DELIM);

            for(String specPosStr : specPositionsStr)
            {
                String[] items = specPosStr.split(":", 2);
                SpecificPositions.add(new BasePosition(items[0], Integer.parseInt(items[1])));
            }

            if(SpecificChrRegions.Regions.isEmpty())
            {
                SpecificPositions.forEach(x -> SpecificChrRegions.addRegion(
                        new ChrBaseRegion(x.Chromosome, x.Position - 300, x.Position + 300)));
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

    public String outputDir() { return pathFromFile(OutputFile); }

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
        if(!makeOutputDir(outputDir))
            return false;

        if(Visualiser.Enabled && !makeOutputDir(new File(Visualiser.OutputDir)))
            return false;

        return true;
    }

    public boolean makeOutputDir(final File outputDir)
    {
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

    public boolean bqrRecordWritingOnly() { return BQR.WriteReads; }

    public boolean logPerfStats() { return PerfWarnTime > 0; }

    public static void registerCommonConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE_BAM, false, REFERENCE_BAMS_DESC);

        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, "Path to sample data files");
        configBuilder.addConfigItem(OUTPUT_VCF, true, "Output VCF filename");

        addRefGenomeConfig(configBuilder, true);

        // is this common?
        configBuilder.addInteger(
                READ_CONTEXT_FLANK_SIZE, "Size of read context flank", DEFAULT_FLANK_LENGTH);

        configBuilder.addInteger(MIN_MAP_QUALITY, "Min map quality to apply to non-hotspot variants", DEFAULT_MIN_MAP_QUALITY);
        configBuilder.addInteger(READ_LENGTH, "Read length, otherwise will sample from BAM", 0);
        configBuilder.addFlag(INCLUDE_MT, "Call MT variants");
        configBuilder.addInteger(SLICE_SIZE, "Slice size", DEFAULT_SLICE_SIZE);
        configBuilder.addInteger(MAX_PARTITION_SLICES, "Max slices per partition", DEFAULT_MAX_PARTITION_SLICES);

        configBuilder.addInteger(MAX_READ_DEPTH, "Max depth to look for evidence", DEFAULT_MAX_READ_DEPTH);
        configBuilder.addInteger(MAX_READ_DEPTH_PANEL, "Max depth to look for evidence in panel", DEFAULT_MAX_READ_DEPTH_PANEL);
        configBuilder.addFlag(NO_FRAGMENT_SYNC, "Disable fragment reads sync in evidence phase");
        configBuilder.addFlag(WRITE_FRAG_LENGTHS, "Write fragment lengths to file");
        addValidationStringencyOption(configBuilder);

        FilterConfig.registerConfig(configBuilder);
        QualityConfig.registerConfig(configBuilder);
        BqrConfig.registerConfig(configBuilder);
        SequencingConfig.registerConfig(configBuilder);
        configBuilder.addPath(JITTER_PARAMS_DIR, false, "Path to sample jitter parameter files");

        VisConfig.registerConfig(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addFlag(LOG_EVIDENCE_READS, "Write evidence read data");

        configBuilder.addFlag(LOG_LPS_DATA, "Log local phasing data");
        configBuilder.addDecimal(PERF_WARN_TIME, "Log details of partitions taking longer than X seconds", 0.0);

        // debug
        configBuilder.addConfigItem(
                SPECIFIC_POSITIONS,
                "Restrict to specific positions(s) of form chromosome:position, separated by ';'");

        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    @VisibleForTesting
    public SageConfig(boolean highDepthMode)
    {
        SampleDataDir = "";
        ReferenceIds = Lists.newArrayList();
        ReferenceBams = Lists.newArrayList();
        Filter = new FilterConfig();
        Quality = new QualityConfig(highDepthMode);
        BQR = new BqrConfig();
        JitterParamsDir = null;
        SpecificChrRegions = new SpecificRegions();
        IncludeMT = false;
        RegionSliceSize = DEFAULT_SLICE_SIZE;
        MinMapQuality = DEFAULT_MIN_MAP_QUALITY;
        MaxReadDepth = DEFAULT_MAX_READ_DEPTH;
        MaxReadDepthPanel = DEFAULT_MAX_READ_DEPTH_PANEL;
        ReadContextFlankLength = DEFAULT_FLANK_LENGTH;
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
        Sequencing = new SequencingConfig(false, SequencingType.ILLUMINA);
        WriteFragmentLengths = false;
        Visualiser = new VisConfig();
        SyncFragments = true;
        SpecificPositions = Collections.emptyList();
        LogEvidenceReads = false;
    }
}
