package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_LEVEL;
import static com.hartwig.hmftools.common.utils.ConfigUtils.containsFlag;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS_DESC;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH_PANEL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MNV;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_SLICE_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.ITEM_DELIM;

import static htsjdk.samtools.ValidationStringency.DEFAULT_STRINGENCY;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.quality.QualityConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.ValidationStringency;

public class SageConfig
{
    public final List<String> ReferenceIds;
    public final List<String> ReferenceBams;
    public final List<String> TumorIds;
    public final List<String> TumorBams;

    public final String SampleDataDir;
    public final String ResourceDir;

    public final boolean AppendMode;
    public final String PanelBed;
    public final boolean PanelOnly;
    public final boolean MnvEnabled;
    public final String Hotspots;
    public final FilterConfig Filter;
    public final QualityConfig Quality;
    public final QualityRecalibrationConfig QualityRecalibration;
    public final Set<String> SpecificChromosomes;
    public final Set<Integer> SpecificPositions;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean IncludeMT;
    public final int RegionSliceSize;
    public final int MinMapQuality;
    public final int MaxReadDepth;
    public final int MaxReadDepthPanel;
    public final int ReadContextFlankSize;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String HighConfidenceBed;
    public final String CoverageBed;
    public final String OutputFile;

    public final String Version;
    public final int Threads;
    public final boolean LogLpsData;
    public final double PerfWarnTime;

    private boolean mIsValid;

    public final ValidationStringency Stringency;

    public int typicalReadLength()
    {
        return 151;
    }

    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String RESOURCE_DIR = "resource_dir";
    
    private static final String THREADS = "threads";
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
    private static final String MNV = "mnv_enabled";
    private static final String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";
    private static final String COVERAGE_BED = "coverage_bed";
    private static final String VALIDATION_STRINGENCY = "validation_stringency";
    private static final String INCLUDE_MT = "include_mt";

    private static final String SPECIFIC_CHROMOSOMES = "specific_chr";
    private static final String SPECIFIC_POSITIONS = "specific_positions";
    private static final String LOG_LPS_DATA = "log_lps_data";
    private static final String PERF_WARN_TIME = "perf_warn_time";

    public SageConfig(boolean appendMode, @NotNull final String version, @NotNull final CommandLine cmd)
    {
        mIsValid = true;
        Version = version;
        
        AppendMode = appendMode;

        RefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

        ReferenceIds = Lists.newArrayList();
        if(cmd.hasOption(REFERENCE))
        {
            ReferenceIds.addAll(Arrays.asList(cmd.getOptionValue(REFERENCE).split(",")));
        }

        TumorIds = Lists.newArrayList();
        if(cmd.hasOption(TUMOR))
        {
            TumorIds.addAll(Arrays.asList(cmd.getOptionValue(TUMOR).split(",")));
        }

        SampleDataDir = cmd.hasOption(SAMPLE_DATA_DIR) ? checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR)) : "";
        ResourceDir = cmd.hasOption(RESOURCE_DIR) ? checkAddDirSeparator(cmd.getOptionValue(RESOURCE_DIR)) : "";

        ReferenceBams = Lists.newArrayList();
        TumorBams = Lists.newArrayList();

        if(cmd.hasOption(REFERENCE_BAM))
        {
            Arrays.stream(cmd.getOptionValue(REFERENCE_BAM, Strings.EMPTY).split(","))
                    .forEach(x -> ReferenceBams.add(SampleDataDir + x));
        }

        if(cmd.hasOption(TUMOR_BAM))
        {
            Arrays.stream(cmd.getOptionValue(TUMOR_BAM, Strings.EMPTY).split(","))
                    .forEach(x -> TumorBams.add(SampleDataDir + x));
        }

        IncludeMT = cmd.hasOption(INCLUDE_MT);

        SpecificChromosomes = Sets.newHashSet();
        SpecificPositions = Sets.newHashSet();
        SpecificRegions = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_REGIONS))
        {
            try
            {
                SpecificRegions.addAll(ChrBaseRegion.loadSpecificRegions(cmd));

                for(ChrBaseRegion region : SpecificRegions)
                {
                    SG_LOGGER.info("filtering for specific region: {}", region);
                    SpecificChromosomes.add(region.Chromosome);
                }
            }
            catch(Exception e)
            {
                SG_LOGGER.error("invalid specific regions: {}", cmd.getOptionValue(SPECIFIC_REGIONS));
                mIsValid = false;
            }
        }
        else if(cmd.hasOption(SPECIFIC_CHROMOSOMES))
        {
            final String chromosomeList = cmd.getOptionValue(SPECIFIC_CHROMOSOMES, Strings.EMPTY);
            if(!chromosomeList.isEmpty())
            {
                SpecificChromosomes.addAll(Lists.newArrayList(chromosomeList.split(ITEM_DELIM)));
            }
        }

        if(cmd.hasOption(SPECIFIC_POSITIONS))
        {
            final String positionList = cmd.getOptionValue(SPECIFIC_POSITIONS, Strings.EMPTY);
            if(!positionList.isEmpty())
            {
                Arrays.stream(positionList.split(ITEM_DELIM)).forEach(x -> SpecificPositions.add(Integer.parseInt(x)));
            }
        }

        OutputFile = SampleDataDir + cmd.getOptionValue(OUTPUT_VCF);

        PanelBed = getReferenceFile(cmd, PANEL_BED);
        CoverageBed = getReferenceFile(cmd, COVERAGE_BED);
        HighConfidenceBed = getReferenceFile(cmd, HIGH_CONFIDENCE_BED);
        Hotspots = getReferenceFile(cmd, HOTSPOTS);
        RefGenomeFile = getReferenceFile(cmd, REF_GENOME);

        Stringency = ValidationStringency.valueOf(cmd.getOptionValue(VALIDATION_STRINGENCY, DEFAULT_STRINGENCY.toString()));
        MnvEnabled = getConfigValue(cmd, MNV, DEFAULT_MNV);
        RegionSliceSize = getConfigValue(cmd, SLICE_SIZE, DEFAULT_SLICE_SIZE);
        ReadContextFlankSize = getConfigValue(cmd, READ_CONTEXT_FLANK_SIZE, DEFAULT_READ_CONTEXT_FLANK_SIZE);
        MinMapQuality = getConfigValue(cmd, MIN_MAP_QUALITY, DEFAULT_MIN_MAP_QUALITY);
        MaxReadDepth = getConfigValue(cmd, MAX_READ_DEPTH, DEFAULT_MAX_READ_DEPTH);
        MaxReadDepthPanel = getConfigValue(cmd, MAX_READ_DEPTH_PANEL, DEFAULT_MAX_READ_DEPTH_PANEL);

        Filter = new FilterConfig(cmd);
        Quality = new QualityConfig(cmd);
        QualityRecalibration = new QualityRecalibrationConfig(cmd);

        PanelOnly = containsFlag(cmd, PANEL_ONLY);
        LogLpsData = containsFlag(cmd, LOG_LPS_DATA);

        PerfWarnTime = Double.parseDouble(cmd.getOptionValue(PERF_WARN_TIME, "0"));

        Threads = getConfigValue(cmd, THREADS, 1);
    }

    private String getReferenceFile(final CommandLine cmd, final String config)
    {
        if(!cmd.hasOption(config))
            return "";

        if(ResourceDir.isEmpty())
            return cmd.getOptionValue(config);

        return ResourceDir + cmd.getOptionValue(config);
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

    public static Options createSageOptions()
    {
        final Options options = new Options();
        options.addOption(MNV, true, "Enable MNVs [" + DEFAULT_MNV + "]");
        options.addOption(TUMOR, true, "Tumor sample, or collection separated by ',\"");
        options.addOption(TUMOR_BAM, true, "Tumor bam file");
        options.addOption(READ_CONTEXT_FLANK_SIZE, true, "Size of read context flank [" + DEFAULT_READ_CONTEXT_FLANK_SIZE + "]");

        options.addOption(HIGH_CONFIDENCE_BED, true, "High confidence regions bed file");
        options.addOption(PANEL_BED, true, "Panel regions bed file");
        options.addOption(PANEL_ONLY, false, "Only examine panel for variants");
        options.addOption(HOTSPOTS, true, "Hotspots");
        options.addOption(COVERAGE_BED, true, "Coverage is calculated for optionally supplied bed");
        options.addOption(VALIDATION_STRINGENCY, true, "SAM validation strategy: STRICT, SILENT, LENIENT [STRICT]");
        options.addOption(LOG_LPS_DATA, false, "Log local phasing data");
        options.addOption(PERF_WARN_TIME, true, "Log details of partitions taking longer than X seconds");

        commonOptions().getOptions().forEach(options::addOption);
        FilterConfig.createOptions().getOptions().forEach(options::addOption);
        addEnsemblDir(options);
        return options;
    }

    public static Options commonOptions()
    {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads [" + 1 + "]");
        options.addOption(REFERENCE, true, "Reference sample, or collection separated by ',");
        options.addOption(RESOURCE_DIR, true, "Resource files");
        options.addOption(SAMPLE_DATA_DIR, true, "Path to sample data files");
        options.addOption(REFERENCE_BAM, true, "Reference bam file");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Assembly, must be one of [37, 38]");
        options.addOption(OUTPUT_VCF, true, "Output vcf");
        options.addOption(MIN_MAP_QUALITY, true, "Min map quality to apply to non-hotspot variants [" + DEFAULT_MIN_MAP_QUALITY + "]");
        options.addOption(SPECIFIC_CHROMOSOMES, true, "Run for subset of chromosomes, split by ';'");
        options.addOption(SPECIFIC_REGIONS, true, SPECIFIC_REGIONS_DESC);
        options.addOption(SPECIFIC_POSITIONS, true, "Run for specific positions(s) separated by ';', for debug purposes");
        options.addOption(INCLUDE_MT, false, "Call MT variants");
        options.addOption(SLICE_SIZE, true, "Slice size [" + DEFAULT_SLICE_SIZE + "]");

        options.addOption(MAX_READ_DEPTH, true, "Max depth to look for evidence [" + DEFAULT_MAX_READ_DEPTH + "]");
        options.addOption(MAX_READ_DEPTH_PANEL, true, "Max depth to look for evidence in panel [" + DEFAULT_MAX_READ_DEPTH_PANEL + "]");

        QualityConfig.createOptions().getOptions().forEach(options::addOption);
        QualityRecalibrationConfig.createOptions().getOptions().forEach(options::addOption);

        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(LOG_LEVEL, true, "Log level");
        return options;
    }

    public String geneCoverageFile(final String sample)
    {
        String filename = sample + ".sage.gene.coverage.tsv";
        String parent = new File(OutputFile).getParent();
        return parent == null ? filename : parent + File.separator + filename;
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
        MnvEnabled = DEFAULT_MNV;
        Hotspots = "hotspots";
        Filter = new FilterConfig();
        Quality = new QualityConfig();
        QualityRecalibration = new QualityRecalibrationConfig();
        SpecificChromosomes = Sets.newHashSet();
        SpecificPositions = Sets.newHashSet();
        SpecificRegions = Lists.newArrayList();
        IncludeMT = false;
        RegionSliceSize = 500_000;
        MinMapQuality = DEFAULT_MIN_MAP_QUALITY;
        MaxReadDepth = DEFAULT_MAX_READ_DEPTH;
        MaxReadDepthPanel = DEFAULT_MAX_READ_DEPTH_PANEL;
        ReadContextFlankSize = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        RefGenomeFile = "refGenome";
        HighConfidenceBed = "highConf";
        CoverageBed = "coverage";
        OutputFile = "out.vcf";
        Version = "1.0";
        Threads = 1;
        LogLpsData = false;
        PerfWarnTime = 0;
        RefGenVersion = V37;
        Stringency = ValidationStringency.DEFAULT_STRINGENCY;
        AppendMode = false;
    }
}
