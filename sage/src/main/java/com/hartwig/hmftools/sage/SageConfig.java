package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
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
import com.hartwig.hmftools.common.samtools.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.quality.QualityConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationConfig;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.ValidationStringency;

public class SageConfig
{
    public final List<String> ReferenceIds;
    public final List<String> ReferenceBams;
    public final List<String> TumorIds;
    public final List<String> TumorBams;

    public final String SampleDataDir;
    public final String ResourceDir;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String HighConfidenceBed;
    public final String CoverageBed;
    public final String OutputFile;
    public final String PanelBed;
    public final boolean PanelOnly;
    public final String Hotspots;
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
    public final int ExpectedReadLength;
    public final ValidationStringency BamStringency;

    public final String Version;
    public final int Threads;

    public final boolean AppendMode;
    public final boolean TrackUMIs;

    // debug
    public final List<String> SpecificChromosomes;
    public final Set<Integer> SpecificPositions;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean LogEvidenceReads;
    public final boolean LogLpsData;
    public final double PerfWarnTime;

    private boolean mIsValid;

    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String RESOURCE_DIR = "resource_dir";
    
    private static final String REFERENCE = "reference";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String TUMOR = "tumor";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String OUTPUT_VCF = "out";
    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String PANEL_BED = "panel_bed";
    private static final String PANEL_ONLY = "panel_only";
    private static final String HOTSPOTS = "hotspots";
    private static final String MAX_READ_DEPTH = "max_read_depth";
    private static final String MAX_READ_DEPTH_PANEL = "max_read_depth_panel";
    private static final String REF_GENOME_VERSION = "ref_genome_version";
    private static final String SLICE_SIZE = "slice_size";
    private static final String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";
    private static final String COVERAGE_BED = "coverage_bed";
    private static final String INCLUDE_MT = "include_mt";
    private static final String EXPECTED_READ_LENGTH = "read_length";
    private static final String SYNC_FRAGMENTS = "sync_fragments";
    private static final String TRACK_UMIS = "track_umis";

    private static final String SPECIFIC_POSITIONS = "specific_positions";
    private static final String LOG_EVIDENCE_READS = "log_evidence_reads";
    private static final String LOG_LPS_DATA = "log_lps_data";
    private static final String PERF_WARN_TIME = "perf_warn_time";

    public SageConfig(boolean appendMode, final String version, final ConfigBuilder configBuilder)
    {
        mIsValid = true;
        Version = version;
        
        AppendMode = appendMode;

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        ReferenceIds = Lists.newArrayList();
        if(configBuilder.hasValue(REFERENCE))
        {
            ReferenceIds.addAll(Arrays.asList(configBuilder.getValue(REFERENCE).split(",")));
        }

        TumorIds = Lists.newArrayList();
        if(!appendMode && configBuilder.hasValue(TUMOR))
        {
            TumorIds.addAll(Arrays.asList(configBuilder.getValue(TUMOR).split(",")));
        }

        SampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR, ""));
        ResourceDir = checkAddDirSeparator(configBuilder.getValue(RESOURCE_DIR, ""));

        ReferenceBams = Lists.newArrayList();
        TumorBams = Lists.newArrayList();

        if(configBuilder.hasValue(REFERENCE_BAM))
        {
            Arrays.stream(configBuilder.getValue(REFERENCE_BAM, Strings.EMPTY).split(","))
                    .forEach(x -> ReferenceBams.add(SampleDataDir + x));
        }

        if(!appendMode && configBuilder.hasValue(TUMOR_BAM))
        {
            Arrays.stream(configBuilder.getValue(TUMOR_BAM, Strings.EMPTY).split(","))
                    .forEach(x -> TumorBams.add(SampleDataDir + x));
        }

        IncludeMT = configBuilder.hasFlag(INCLUDE_MT);

        SpecificPositions = Sets.newHashSet();
        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        try
        {
            loadSpecificChromsomesOrRegions(configBuilder, SpecificChromosomes, SpecificRegions, SG_LOGGER);
        }
        catch(ParseException e)
        {
            mIsValid = false;
        }

        if(configBuilder.hasValue(SPECIFIC_POSITIONS))
        {
            final String positionList = configBuilder.getValue(SPECIFIC_POSITIONS, Strings.EMPTY);
            if(!positionList.isEmpty())
            {
                Arrays.stream(positionList.split(ITEM_DELIM)).forEach(x -> SpecificPositions.add(Integer.parseInt(x)));
            }
        }

        OutputFile = SampleDataDir + configBuilder.getValue(OUTPUT_VCF);

        PanelBed = getReferenceFile(configBuilder, PANEL_BED);
        CoverageBed = getReferenceFile(configBuilder, COVERAGE_BED);
        HighConfidenceBed = getReferenceFile(configBuilder, HIGH_CONFIDENCE_BED);
        Hotspots = getReferenceFile(configBuilder, HOTSPOTS);
        RefGenomeFile = getReferenceFile(configBuilder, REF_GENOME);

        BamStringency = BamUtils.validationStringency(configBuilder);
        RegionSliceSize = configBuilder.getInteger(SLICE_SIZE);
        ReadContextFlankSize = configBuilder.getInteger(READ_CONTEXT_FLANK_SIZE);
        MinMapQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxReadDepth = configBuilder.getInteger(MAX_READ_DEPTH);
        MaxReadDepthPanel = configBuilder.getInteger(MAX_READ_DEPTH_PANEL);
        ExpectedReadLength = configBuilder.getInteger(EXPECTED_READ_LENGTH);
        SyncFragments = configBuilder.hasFlag(SYNC_FRAGMENTS);

        Filter = new FilterConfig(configBuilder);
        Quality = new QualityConfig(configBuilder);
        QualityRecalibration = new QualityRecalibrationConfig(configBuilder);

        TrackUMIs = configBuilder.hasFlag(TRACK_UMIS);

        PanelOnly = configBuilder.hasFlag(PANEL_ONLY);
        LogLpsData = configBuilder.hasFlag(LOG_LPS_DATA);
        LogEvidenceReads = !SpecificRegions.isEmpty() && configBuilder.hasFlag(LOG_EVIDENCE_READS);

        if(LogEvidenceReads)
        {
            SG_LOGGER.trace("READ_EV,SampleId,Chromosome,Position,Ref,Alt,MatchType,ReadId,ReadStart,Cigar,LeftCore,Index,RightCore,ReadIndex");
        }

        PerfWarnTime = configBuilder.getDecimal(PERF_WARN_TIME);

        Threads = parseThreads(configBuilder);
    }

    private String getReferenceFile(final ConfigBuilder configBuilder, final String config)
    {
        if(!configBuilder.hasValue(config))
            return "";

        if(ResourceDir.isEmpty())
            return configBuilder.getValue(config);

        return ResourceDir + configBuilder.getValue(config);
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

        if(TumorIds.size() != TumorBams.size())
        {
            SG_LOGGER.error("Each tumor sample must have matching bam");
            return false;
        }

        for(String tumorBam : TumorBams)
        {
            if(!new File(tumorBam).exists())
            {
                SG_LOGGER.error("Unable to locate tumor bam({})", tumorBam);
                return false;
            }
        }

        final File outputDir = new File(OutputFile).getParentFile();
        if(outputDir != null && !outputDir.exists() && !outputDir.mkdirs())
        {
            SG_LOGGER.error("unable to write directory({})", outputDir.toString());
            return false;
        }

        if(!AppendMode)
        {
            if(TumorIds.isEmpty())
            {
                SG_LOGGER.error("At least one tumor must be supplied");
                return false;
            }
        }

        return true;
    }

    public boolean processChromosome(final String chromosome)
    {
        if(!SpecificChromosomes.isEmpty() && !SpecificChromosomes.contains(chromosome))
            return false;

        if(HumanChromosome.contains(chromosome))
            return true;

        if(IncludeMT && MitochondrialChromosome.contains(chromosome))
            return true;

        return false;
    }

    public boolean logPerfStats() { return PerfWarnTime > 0; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, true, "Tumor sample, or collection separated by ','");
        configBuilder.addConfigItem(TUMOR_BAM, true, "Tumor bam file(s)");

        configBuilder.addInteger(
                READ_CONTEXT_FLANK_SIZE, "Size of read context flank", DEFAULT_READ_CONTEXT_FLANK_SIZE);

        configBuilder.addConfigItem(HIGH_CONFIDENCE_BED, false, "High confidence regions bed file");
        configBuilder.addConfigItem(PANEL_BED, false, "Panel regions bed file");
        configBuilder.addConfigItem(PANEL_ONLY, false, "Only examine panel for variants");
        configBuilder.addConfigItem(HOTSPOTS, false, "Hotspots");
        configBuilder.addConfigItem(COVERAGE_BED, false, "Coverage is calculated for optionally supplied bed");
        configBuilder.addFlag(LOG_LPS_DATA, "Log local phasing data");
        configBuilder.addDecimal(PERF_WARN_TIME, "Log details of partitions taking longer than X seconds", 0.0);

        registerCommonConfig(configBuilder);
        FilterConfig.registerConfig(configBuilder);
        addEnsemblDir(configBuilder);
    }

    public static void registerCommonConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REFERENCE, false, "Reference sample, or collection separated by ','");
        configBuilder.addConfigItem(REFERENCE_BAM, false, "Reference bam file");

        configBuilder.addPath(RESOURCE_DIR, false, "Resource files");
        configBuilder.addPath(SAMPLE_DATA_DIR, false, "Path to sample data files");
        configBuilder.addConfigItem(OUTPUT_VCF, true, "Output vcf");

        // addRefGenomeConfig(configBuilder, true);
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        configBuilder.addConfigItem(REF_GENOME, true, REF_GENOME_CFG_DESC);


        configBuilder.addInteger(MIN_MAP_QUALITY, "Min map quality to apply to non-hotspot variants", DEFAULT_MIN_MAP_QUALITY);
        configBuilder.addInteger(EXPECTED_READ_LENGTH, "Expected read length", DEFAULT_READ_LENGTH);
        configBuilder.addFlag(INCLUDE_MT, "Call MT variants");
        configBuilder.addInteger(SLICE_SIZE, "Slice size", DEFAULT_SLICE_SIZE);

        configBuilder.addInteger(MAX_READ_DEPTH, "Max depth to look for evidence", DEFAULT_MAX_READ_DEPTH);
        configBuilder.addInteger(MAX_READ_DEPTH_PANEL, "Max depth to look for evidence in panel", DEFAULT_MAX_READ_DEPTH_PANEL);
        configBuilder.addFlag(SYNC_FRAGMENTS, "Handle overlapping fragment reads in evidence phase");
        configBuilder.addFlag(TRACK_UMIS, "Record counts of UMI types");
        addValidationStringencyOption(configBuilder);

        QualityConfig.registerConfig(configBuilder);
        QualityRecalibrationConfig.registerConfig(configBuilder);

        // debug
        configBuilder.addConfigItem(SPECIFIC_POSITIONS, "Run for specific positions(s) separated by ';', for debug purposes");
        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addFlag(LOG_EVIDENCE_READS, "Write evidence read data");

        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    public SageConfig()
    {
        SampleDataDir = "";
        ResourceDir = "";
        ReferenceIds = Lists.newArrayList();
        ReferenceBams = Lists.newArrayList();
        TumorIds = Lists.newArrayList();
        TumorBams = Lists.newArrayList();
        PanelBed = "panel";
        PanelOnly = false;
        Hotspots = "hotspots";
        Filter = new FilterConfig();
        Quality = new QualityConfig();
        QualityRecalibration = new QualityRecalibrationConfig();
        SpecificChromosomes = Lists.newArrayList();
        SpecificPositions = Sets.newHashSet();
        SpecificRegions = Lists.newArrayList();
        IncludeMT = false;
        RegionSliceSize = 500_000;
        MinMapQuality = DEFAULT_MIN_MAP_QUALITY;
        MaxReadDepth = DEFAULT_MAX_READ_DEPTH;
        MaxReadDepthPanel = DEFAULT_MAX_READ_DEPTH_PANEL;
        ReadContextFlankSize = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        ExpectedReadLength = DEFAULT_READ_LENGTH;
        RefGenomeFile = "refGenome";
        HighConfidenceBed = "highConf";
        CoverageBed = "coverage";
        OutputFile = "out.vcf";
        Version = "1.0";
        Threads = 1;
        LogLpsData = false;
        PerfWarnTime = 0;
        RefGenVersion = V37;
        BamStringency = ValidationStringency.DEFAULT_STRINGENCY;
        AppendMode = false;
        TrackUMIs = false;
        SyncFragments = false;
        LogEvidenceReads = false;
    }
}
