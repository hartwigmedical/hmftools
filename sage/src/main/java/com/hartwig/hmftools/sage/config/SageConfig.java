package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.containsFlag;
import static com.hartwig.hmftools.common.utils.ConfigUtils.defaultEnumValue;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

import java.io.File;
import java.io.IOException;
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
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.ValidationStringency;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SageConfig
{
    String THREADS = "threads";
    String REFERENCE = "reference";
    String REFERENCE_BAM = "reference_bam";
    String TUMOR = "tumor";
    String TUMOR_BAM = "tumor_bam";
    String REF_GENOME = "ref_genome";
    String OUTPUT_VCF = "out";
    String INPUT_VCF = "input_vcf";
    String MIN_MAP_QUALITY = "min_map_quality";
    String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    String PANEL_BED = "panel_bed";
    String PANEL_ONLY = "panel_only";
    String HOTSPOTS = "hotspots";
    String MAX_READ_DEPTH = "max_read_depth";
    String MAX_READ_DEPTH_PANEL = "max_read_depth_panel";
    String MAX_REALIGNMENT_DEPTH = "max_realignment_depth";
    String REF_GENOME_VERSION = "ref_genome_version";
    String CHR = "chr";
    String SLICE_SIZE = "slice_size";
    String MNV = "mnv_enabled";
    String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";
    String COVERAGE_BED = "coverage_bed";
    String VALIDATION_STRINGENCY = "validation_stringency";

    int DEFAULT_THREADS = 2;
    int DEFAULT_MIN_MAP_QUALITY = 10;
    int DEFAULT_MAX_READ_DEPTH = 1000;
    int DEFAULT_MAX_READ_DEPTH_PANEL = 100_000;
    int DEFAULT_MAX_REALIGNMENT_DEPTH = 1000;
    int DEFAULT_SLICE_SIZE = 100_000;
    int DEFAULT_READ_CONTEXT_FLANK_SIZE = 10;
    boolean DEFAULT_MNV = true;

    @NotNull
    static Options createSageOptions()
    {
        final Options options = new Options();
        options.addOption(MNV, true, "Enable MNVs [" + DEFAULT_MNV + "]");
        options.addOption(REF_GENOME_VERSION, true, "Assembly, must be one of [37, 38]");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
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

    @NotNull
    static Options createAddReferenceOptions()
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
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(REF_GENOME, true, "Path to indexed ref genome fasta file");
        options.addOption(OUTPUT_VCF, true, "Path to output vcf");
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

    @NotNull
    String version();

    int threads();

    @NotNull
    List<String> reference();

    @NotNull
    List<String> referenceBam();

    @NotNull
    String refGenome();

    @NotNull
    String highConfidenceBed();

    @NotNull
    String coverageBed();

    @NotNull
    List<String> tumor();

    @NotNull
    List<String> tumorBam();

    @NotNull
    String inputFile();

    @NotNull
    String outputFile();

    @NotNull
    List<HmfTranscriptRegion> transcriptRegions();

    @NotNull
    ValidationStringency validationStringency();

    @NotNull
    default String geneCoverageFile(@NotNull final String sample)
    {
        String filename = sample + ".sage.gene.coverage.tsv";
        String parent = new File(outputFile()).getParent();
        return parent == null ? filename : parent + File.separator + filename;
    }

    @NotNull
    default String baseQualityRecalibrationFile(@NotNull final String sample)
    {
        String parent = new File(outputFile()).getParent();
        return parent == null ? sample + ".sage.bqr.tsv" : parent + File.separator + sample + ".sage.bqr.tsv";
    }

    @NotNull
    String panelBed();

    boolean panelOnly();

    boolean mnvEnabled();

    @NotNull
    String hotspots();

    @NotNull
    FilterConfig filter();

    @NotNull
    QualityConfig qualityConfig();

    @NotNull
    BaseQualityRecalibrationConfig baseQualityRecalibrationConfig();

    @NotNull
    Set<String> chromosomes();

    default int typicalReadLength()
    {
        return 151;
    }

    int regionSliceSize();

    int minMapQuality();

    int maxRealignmentDepth();

    int maxReadDepth();

    int maxReadDepthPanel();

    default int maxSkippedReferenceRegions()
    {
        return 50;
    }

    int readContextFlankSize();

    @NotNull
    static SageConfig createConfig(boolean appendMode, @NotNull final String version, @NotNull final CommandLine cmd)
            throws ParseException, IOException
    {
        final int threads = getConfigValue(cmd, THREADS, DEFAULT_THREADS);
        final RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION));

        final List<String> referenceList = Lists.newArrayList();
        if(cmd.hasOption(REFERENCE))
        {
            referenceList.addAll(Arrays.asList(cmd.getOptionValue(REFERENCE).split(",")));
        }

        final List<String> referenceBamList = Lists.newArrayList();
        if(cmd.hasOption(REFERENCE_BAM))
        {
            referenceBamList.addAll(Arrays.asList(cmd.getOptionValue(REFERENCE_BAM, Strings.EMPTY).split(",")));
        }

        if(referenceList.size() != referenceBamList.size())
        {
            throw new ParseException("Each reference sample must have matching bam");
        }

        for(String referenceBam : referenceBamList)
        {
            if(!new File(referenceBam).exists())
            {
                throw new ParseException("Unable to locate reference bam " + referenceBam);
            }
        }

        if(!cmd.hasOption(REF_GENOME))
        {
            throw new ParseException(REF_GENOME + " is a mandatory argument");
        }

        if(!cmd.hasOption(OUTPUT_VCF))
        {
            throw new ParseException(OUTPUT_VCF + " is a mandatory argument");
        }

        final List<String> tumorList = Lists.newArrayList();
        if(cmd.hasOption(TUMOR))
        {
            tumorList.addAll(Arrays.asList(cmd.getOptionValue(TUMOR).split(",")));
        }

        final List<String> tumorBamList = Lists.newArrayList();
        if(cmd.hasOption(TUMOR_BAM))
        {
            tumorBamList.addAll(Arrays.asList(cmd.getOptionValue(TUMOR_BAM, Strings.EMPTY).split(",")));
        }

        if(tumorList.size() != tumorBamList.size())
        {
            throw new ParseException("Each tumor sample must have matching bam");
        }

        for(String tumorBam : tumorBamList)
        {
            if(!new File(tumorBam).exists())
            {
                throw new ParseException("Unable to locate tumor bam " + tumorBam);
            }
        }

        final Set<String> chromosomes = Sets.newHashSet();
        final String chromosomeList = cmd.getOptionValue(CHR, Strings.EMPTY);
        if(!chromosomeList.isEmpty())
        {
            chromosomes.addAll(Lists.newArrayList(chromosomeList.split(",")));
        }

        final List<HmfTranscriptRegion> transcripts;
        final String outputVcf = cmd.getOptionValue(OUTPUT_VCF);
        final String inputVcf = cmd.getOptionValue(INPUT_VCF, Strings.EMPTY);
        if(appendMode)
        {
            transcripts = Lists.newArrayList();
            if(!cmd.hasOption(INPUT_VCF))
            {
                throw new ParseException(INPUT_VCF + " is a mandatory argument");
            }

            if(outputVcf.equals(inputVcf))
            {
                throw new ParseException("Input and output VCFs must be different");
            }

            if(referenceList.isEmpty())
            {
                throw new ParseException("At least one reference must be supplied");
            }

        }
        else
        {

            if(tumorList.isEmpty())
            {
                throw new ParseException("At least one tumor must be supplied");
            }

            transcripts = refGenomeVersion.is37() ? HmfGenePanelSupplier.allGeneList37() : HmfGenePanelSupplier.allGeneList38();
        }

        final File outputDir = new File(outputVcf).getParentFile();
        if(outputDir != null && !outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write directory " + outputDir.toString());
        }

        final ValidationStringency validationStringency =
                defaultEnumValue(cmd, VALIDATION_STRINGENCY, ValidationStringency.DEFAULT_STRINGENCY);

        return ImmutableSageConfig.builder()
                .version(version)
                .transcriptRegions(transcripts)
                .inputFile(inputVcf)
                .outputFile(outputVcf)
                .threads(threads)
                .chromosomes(chromosomes)
                .reference(referenceList)
                .referenceBam(referenceBamList)
                .tumor(tumorList)
                .tumorBam(tumorBamList)
                .mnvEnabled(getConfigValue(cmd, MNV, DEFAULT_MNV))
                .refGenome(cmd.getOptionValue(REF_GENOME))
                .regionSliceSize(getConfigValue(cmd, SLICE_SIZE, DEFAULT_SLICE_SIZE))
                .readContextFlankSize(getConfigValue(cmd, READ_CONTEXT_FLANK_SIZE, DEFAULT_READ_CONTEXT_FLANK_SIZE))
                .minMapQuality(getConfigValue(cmd, MIN_MAP_QUALITY, DEFAULT_MIN_MAP_QUALITY))
                .maxReadDepth(getConfigValue(cmd, MAX_READ_DEPTH, DEFAULT_MAX_READ_DEPTH))
                .maxReadDepthPanel(getConfigValue(cmd, MAX_READ_DEPTH_PANEL, DEFAULT_MAX_READ_DEPTH_PANEL))
                .maxRealignmentDepth(getConfigValue(cmd, MAX_REALIGNMENT_DEPTH, DEFAULT_MAX_REALIGNMENT_DEPTH))
                .filter(FilterConfig.createConfig(cmd))
                .panelBed(cmd.getOptionValue(PANEL_BED, Strings.EMPTY))
                .coverageBed(cmd.getOptionValue(COVERAGE_BED, Strings.EMPTY))
                .highConfidenceBed(cmd.getOptionValue(HIGH_CONFIDENCE_BED, Strings.EMPTY))
                .hotspots(cmd.getOptionValue(HOTSPOTS, Strings.EMPTY))
                .qualityConfig(QualityConfig.createConfig(cmd, transcripts))
                .baseQualityRecalibrationConfig(BaseQualityRecalibrationConfig.createConfig(cmd))
                .panelOnly(containsFlag(cmd, PANEL_ONLY))
                .validationStringency(validationStringency)
                .build();
    }
}
