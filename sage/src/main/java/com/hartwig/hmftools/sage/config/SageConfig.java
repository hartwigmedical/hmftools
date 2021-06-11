package com.hartwig.hmftools.sage.config;

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

import static htsjdk.samtools.ValidationStringency.DEFAULT_STRINGENCY;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

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

    public final String PanelBed;
    public final boolean PanelOnly;
    public final boolean MnvEnabled;
    public final String Hotspots;
    public final FilterConfig Filter;
    public final QualityConfig Quality;
    public final BaseQualityRecalibrationConfig BaseQualityRecalibration;
    public final Set<String> Chromosomes;
    public final int RegionSliceSize;
    public final int MinMapQuality;
    public final int MaxRealignmentDepth;
    public final int MaxReadDepth;
    public final int MaxReadDepthPanel;
    public final int ReadContextFlankSize;

    public final String RefGenomeFile;
    public final String HighConfidenceBed;
    public final String CoverageBed;
    public final String InputFile;
    public final String OutputFile;

    public final String Version;
    public final int Threads;

    public final List<HmfTranscriptRegion> TranscriptRegions;

    public final ValidationStringency Stringency;

    private final boolean mAppendMode;

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
    private static final String REF_GENOME = "ref_genome";
    private static final String OUTPUT_VCF = "out";
    private static final String INPUT_VCF = "input_vcf";
    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String PANEL_BED = "panel_bed";
    private static final String PANEL_ONLY = "panel_only";
    private static final String HOTSPOTS = "hotspots";
    private static final String MAX_READ_DEPTH = "max_read_depth";
    private static final String MAX_READ_DEPTH_PANEL = "max_read_depth_panel";
    private static final String MAX_REALIGNMENT_DEPTH = "max_realignment_depth";
    private static final String REF_GENOME_VERSION = "ref_genome_version";
    private static final String CHR = "chr";
    private static final String SLICE_SIZE = "slice_size";
    private static final String MNV = "mnv_enabled";
    private static final String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";
    private static final String COVERAGE_BED = "coverage_bed";
    private static final String VALIDATION_STRINGENCY = "validation_stringency";

    public SageConfig(boolean appendMode, @NotNull final String version, @NotNull final CommandLine cmd)
    {
        Version = version;
        
        mAppendMode = appendMode;

        final RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION));

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
        final String chromosomeList = cmd.getOptionValue(CHR, Strings.EMPTY);
        if(!chromosomeList.isEmpty())
        {
            Chromosomes.addAll(Lists.newArrayList(chromosomeList.split(",")));
        }

        TranscriptRegions = Lists.newArrayList();

        if(!appendMode)
        {
            TranscriptRegions.addAll(refGenomeVersion.is37() ? HmfGenePanelSupplier.allGeneList37() : HmfGenePanelSupplier.allGeneList38());
        }

        OutputFile = SampleDataDir + cmd.getOptionValue(OUTPUT_VCF);
        InputFile = SampleDataDir + cmd.getOptionValue(INPUT_VCF, "");

        PanelBed = ResourceDir + cmd.getOptionValue(PANEL_BED, Strings.EMPTY);
        CoverageBed = ResourceDir + cmd.getOptionValue(COVERAGE_BED, Strings.EMPTY);
        HighConfidenceBed = ResourceDir + cmd.getOptionValue(HIGH_CONFIDENCE_BED, Strings.EMPTY);
        Hotspots = ResourceDir +  cmd.getOptionValue(HOTSPOTS, Strings.EMPTY);
        RefGenomeFile = ResourceDir + cmd.getOptionValue(REF_GENOME);

        Stringency = ValidationStringency.valueOf(cmd.getOptionValue(VALIDATION_STRINGENCY, DEFAULT_STRINGENCY.toString()));
        MnvEnabled = getConfigValue(cmd, MNV, DEFAULT_MNV);
        RegionSliceSize = getConfigValue(cmd, SLICE_SIZE, DEFAULT_SLICE_SIZE);
        ReadContextFlankSize = getConfigValue(cmd, READ_CONTEXT_FLANK_SIZE, DEFAULT_READ_CONTEXT_FLANK_SIZE);
        MinMapQuality = getConfigValue(cmd, MIN_MAP_QUALITY, DEFAULT_MIN_MAP_QUALITY);
        MaxReadDepth = getConfigValue(cmd, MAX_READ_DEPTH, DEFAULT_MAX_READ_DEPTH);
        MaxReadDepthPanel = getConfigValue(cmd, MAX_READ_DEPTH_PANEL, DEFAULT_MAX_READ_DEPTH_PANEL);
        MaxRealignmentDepth = getConfigValue(cmd, MAX_REALIGNMENT_DEPTH, DEFAULT_MAX_REALIGNMENT_DEPTH);

        Filter = new FilterConfig(cmd);
        Quality = new QualityConfig(cmd, TranscriptRegions);
        BaseQualityRecalibration = new BaseQualityRecalibrationConfig(cmd);

        PanelOnly = containsFlag(cmd, PANEL_ONLY);

        Threads = getConfigValue(cmd, THREADS, DEFAULT_THREADS);
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
            SG_LOGGER.error("Unable to write directory({})", outputDir.toString());
            return false;
        }

        if(mAppendMode)
        {
            if(InputFile.isEmpty())
            {
                SG_LOGGER.error("No input VCF file specified");
                return false;
            }

            if(InputFile.equals(OutputFile))
            {
                SG_LOGGER.error("Input and output VCFs must be different");
                return false;
            }

            if(ReferenceIds.isEmpty())
            {
                SG_LOGGER.error("At least one reference must be supplied");
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
        options.addOption(REF_GENOME, true, "Indexed ref genome fasta file");
        options.addOption(OUTPUT_VCF, true, "Output vcf");
        options.addOption(MIN_MAP_QUALITY, true, "Min map quality to apply to non-hotspot variants [" + DEFAULT_MIN_MAP_QUALITY + "]");
        options.addOption(CHR, true, "Run for single chromosome");
        options.addOption(SLICE_SIZE, true, "Slice size [" + DEFAULT_SLICE_SIZE + "]");

        options.addOption(MAX_READ_DEPTH, true, "Max depth to look for evidence [" + DEFAULT_MAX_READ_DEPTH + "]");
        options.addOption(MAX_READ_DEPTH_PANEL, true, "Max depth to look for evidence in panel [" + DEFAULT_MAX_READ_DEPTH_PANEL + "]");
        options.addOption(MAX_REALIGNMENT_DEPTH, true, "Max depth to check for realignment [" + DEFAULT_MAX_REALIGNMENT_DEPTH + "]");
        QualityConfig.createOptions().getOptions().forEach(options::addOption);
        BaseQualityRecalibrationConfig.createOptions().getOptions().forEach(options::addOption);
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
        BaseQualityRecalibration = new BaseQualityRecalibrationConfig();
        Chromosomes = Sets.newHashSet();
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
        Version = "1.0";
        Threads = DEFAULT_THREADS;
        TranscriptRegions = HmfGenePanelSupplier.allGeneList37();
        Stringency = ValidationStringency.DEFAULT_STRINGENCY;
        mAppendMode = false;
    }
}
