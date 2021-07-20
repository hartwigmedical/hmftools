package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_LEVEL;
import static com.hartwig.hmftools.common.utils.ConfigUtils.containsFlag;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH_PANEL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_REALIGNMENT_DEPTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MNV;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_SLICE_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_THREADS;
import static com.hartwig.hmftools.sage.SageConstants.ITEM_DELIM;
import static com.hartwig.hmftools.sage.SageConstants.SUB_ITEM_DELIM;

import static htsjdk.samtools.ValidationStringency.DEFAULT_STRINGENCY;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
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
    public final Set<String> Chromosomes;
    public final List<BaseRegion> SpecificRegions;
    public final int RegionSliceSize;
    public final int MinMapQuality;
    public final int MaxRealignmentDepth;
    public final int MaxReadDepth;
    public final int MaxReadDepthPanel;
    public final int ReadContextFlankSize;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String HighConfidenceBed;
    public final String CoverageBed;
    public final String InputFile;
    public final String OutputFile;
    public final boolean WriteCsv;

    public final String Version;
    public final int Threads;

    public final ValidationStringency Stringency;

    public int typicalReadLength()
    {
        return 151;
    }
    public int maxSkippedReferenceRegions() { return 50; }

    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String RESOURCE_DIR = "resource_dir";
    
    private static final String THREADS = "threads";
    private static final String REFERENCE = "reference";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String TUMOR = "tumor";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String OUTPUT_VCF = "out";
    private static final String INPUT_VCF = "input_vcf";
    private static final String WRITE_CSV = "write_csv";
    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String PANEL_BED = "panel_bed";
    private static final String PANEL_ONLY = "panel_only";
    private static final String HOTSPOTS = "hotspots";
    private static final String MAX_READ_DEPTH = "max_read_depth";
    private static final String MAX_READ_DEPTH_PANEL = "max_read_depth_panel";
    private static final String MAX_REALIGNMENT_DEPTH = "max_realignment_depth";
    private static final String REF_GENOME_VERSION = "ref_genome_version";
    private static final String SPECIFIC_CHROMOSOMES = "specific_chr";
    private static final String SPECIFIC_REGIONS = "specific_regions";
    private static final String SLICE_SIZE = "slice_size";
    private static final String MNV = "mnv_enabled";
    private static final String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";
    private static final String COVERAGE_BED = "coverage_bed";
    private static final String VALIDATION_STRINGENCY = "validation_stringency";

    public SageConfig(boolean appendMode, @NotNull final String version, @NotNull final CommandLine cmd)
    {
        Version = version;
        
        AppendMode = appendMode;

        RefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION));

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

        Chromosomes = Sets.newHashSet();
        SpecificRegions = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_REGIONS))
        {
            final List<String> regionStrs = Arrays.stream(cmd.getOptionValue(SPECIFIC_REGIONS).split(ITEM_DELIM, -1)).collect(Collectors.toList());
            for(String regionStr : regionStrs)
            {
                final String[] items = regionStr.split(SUB_ITEM_DELIM);
                if(items.length == 3)
                {
                    BaseRegion region = new BaseRegion(items[0], Integer.parseInt(items[1]), Integer.parseInt(items[2]));

                    if(!region.isValid())
                    {
                        SG_LOGGER.error("invalid specific region: {}", region);
                        continue;
                    }

                    SG_LOGGER.info("filtering for specific region: {}", region);
                    SpecificRegions.add(region);
                    Chromosomes.add(region.Chromosome);
                }
            }
        }
        else if(cmd.hasOption(SPECIFIC_CHROMOSOMES))
        {
            final String chromosomeList = cmd.getOptionValue(SPECIFIC_CHROMOSOMES, Strings.EMPTY);
            if(!chromosomeList.isEmpty())
            {
                Chromosomes.addAll(Lists.newArrayList(chromosomeList.split(ITEM_DELIM)));
            }
        }

        OutputFile = SampleDataDir + cmd.getOptionValue(OUTPUT_VCF);
        WriteCsv = cmd.hasOption(WRITE_CSV);

        InputFile = SampleDataDir + cmd.getOptionValue(INPUT_VCF, "");

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
        MaxRealignmentDepth = getConfigValue(cmd, MAX_REALIGNMENT_DEPTH, DEFAULT_MAX_REALIGNMENT_DEPTH);

        Filter = new FilterConfig(cmd);
        Quality = new QualityConfig(cmd);
        QualityRecalibration = new QualityRecalibrationConfig(cmd);

        PanelOnly = containsFlag(cmd, PANEL_ONLY);

        Threads = getConfigValue(cmd, THREADS, DEFAULT_THREADS);
    }

    private String getReferenceFile(final CommandLine cmd, final String config)
    {
        if(!cmd.hasOption(config))
            return "";

        if(ResourceDir.isEmpty())
            return cmd.getOptionValue(config);

        return ResourceDir + cmd.getOptionValue(config);
    }

    public boolean isValid()
    {
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

        if(AppendMode)
        {
            if(InputFile.isEmpty())
            {
                SG_LOGGER.error("no input VCF file specified");
                return false;
            }

            if(InputFile.equals(OutputFile))
            {
                SG_LOGGER.error("input and output VCFs must be different");
                return false;
            }

            if(ReferenceIds.isEmpty())
            {
                SG_LOGGER.error("at least one reference must be supplied");
                return false;
            }
        }
        else
        {

            if(TumorIds.isEmpty())
            {
                SG_LOGGER.error("At least one tumor must be supplied");
                return false;
            }
        }

        return true;
    }

    public static Options createSageOptions()
    {
        final Options options = new Options();
        options.addOption(MNV, true, "Enable MNVs [" + DEFAULT_MNV + "]");
        options.addOption(REF_GENOME_VERSION, true, "Assembly, must be one of [37, 38]");
        options.addOption(TUMOR, true, "Tumor sample, or collection separated by ',\"");
        options.addOption(TUMOR_BAM, true, "Tumor bam file");
        options.addOption(READ_CONTEXT_FLANK_SIZE, true, "Size of read context flank [" + DEFAULT_READ_CONTEXT_FLANK_SIZE + "]");

        options.addOption(HIGH_CONFIDENCE_BED, true, "High confidence regions bed file");
        options.addOption(PANEL_BED, true, "Panel regions bed file");
        options.addOption(PANEL_ONLY, false, "Only examine panel for variants");
        options.addOption(HOTSPOTS, true, "Hotspots");
        options.addOption(COVERAGE_BED, true, "Coverage is calculated for optionally supplied bed");
        options.addOption(VALIDATION_STRINGENCY, true, "SAM validation strategy: STRICT, SILENT, LENIENT [STRICT]");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(LOG_LEVEL, true, "Log level");

        commonOptions().getOptions().forEach(options::addOption);
        FilterConfig.createOptions().getOptions().forEach(options::addOption);
        return options;
    }

    public static Options createAddReferenceOptions()
    {
        final Options options = new Options();
        options.addOption(INPUT_VCF, true, "Path to input vcf");
        commonOptions().getOptions().forEach(options::addOption);
        return options;
    }

    static Options commonOptions()
    {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Reference sample, or collection separated by ',");
        options.addOption(RESOURCE_DIR, true, "Resource files");
        options.addOption(SAMPLE_DATA_DIR, true, "Path to sample data files");
        options.addOption(REFERENCE_BAM, true, "Reference bam file");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(OUTPUT_VCF, true, "Output vcf");
        options.addOption(WRITE_CSV, false, "Write variant data to CSV as well as VCF");
        options.addOption(MIN_MAP_QUALITY, true, "Min map quality to apply to non-hotspot variants [" + DEFAULT_MIN_MAP_QUALITY + "]");
        options.addOption(SPECIFIC_CHROMOSOMES, true, "Run for subset of chromosomes, split by ';'");
        options.addOption(SPECIFIC_REGIONS, true, "Run for specific regions(s) separated by ';' in format Chr:PosStart:PosEnd");
        options.addOption(SLICE_SIZE, true, "Slice size [" + DEFAULT_SLICE_SIZE + "]");

        options.addOption(MAX_READ_DEPTH, true, "Max depth to look for evidence [" + DEFAULT_MAX_READ_DEPTH + "]");
        options.addOption(MAX_READ_DEPTH_PANEL, true, "Max depth to look for evidence in panel [" + DEFAULT_MAX_READ_DEPTH_PANEL + "]");
        options.addOption(MAX_REALIGNMENT_DEPTH, true, "Max depth to check for realignment [" + DEFAULT_MAX_REALIGNMENT_DEPTH + "]");
        QualityConfig.createOptions().getOptions().forEach(options::addOption);
        QualityRecalibrationConfig.createOptions().getOptions().forEach(options::addOption);
        return options;
    }

    public String geneCoverageFile(@NotNull final String sample)
    {
        String filename = sample + ".sage.gene.coverage.tsv";
        String parent = new File(OutputFile).getParent();
        return parent == null ? filename : parent + File.separator + filename;
    }

    public String baseQualityRecalibrationFile(@NotNull final String sample)
    {
        String parent = new File(OutputFile).getParent();
        return parent == null ? sample + ".sage.bqr.tsv" : parent + File.separator + sample + ".sage.bqr.tsv";
    }

    public SageConfig()
    {
        SampleDataDir = "";
        ResourceDir = "";
        ReferenceIds = Lists.newArrayList("referencIds");
        ReferenceBams = Lists.newArrayList("referenceBams");
        TumorIds = Lists.newArrayList("tumorIds");
        TumorBams = Lists.newArrayList("tumorBams");
        PanelBed = "panel";
        PanelOnly = false;
        MnvEnabled = DEFAULT_MNV;
        Hotspots = "hotspots";
        Filter = new FilterConfig();
        Quality = new QualityConfig();
        QualityRecalibration = new QualityRecalibrationConfig();
        Chromosomes = Sets.newHashSet();
        SpecificRegions = Lists.newArrayList();
        RegionSliceSize = 500_000;
        MinMapQuality = DEFAULT_MIN_MAP_QUALITY;
        MaxRealignmentDepth = DEFAULT_MAX_REALIGNMENT_DEPTH;
        MaxReadDepth = DEFAULT_MAX_READ_DEPTH;
        MaxReadDepthPanel = DEFAULT_MAX_READ_DEPTH_PANEL;
        ReadContextFlankSize = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        RefGenomeFile = "refGenome";
        HighConfidenceBed = "highConf";
        CoverageBed = "coverage";
        InputFile = "in.vcf";
        OutputFile = "out.vcf";
        WriteCsv = false;
        Version = "1.0";
        Threads = DEFAULT_THREADS;
        RefGenVersion = V37;
        Stringency = ValidationStringency.DEFAULT_STRINGENCY;
        AppendMode = false;
    }
}
