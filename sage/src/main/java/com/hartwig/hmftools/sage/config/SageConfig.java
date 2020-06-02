package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.cli.Configs.defaultBooleanValue;
import static com.hartwig.hmftools.common.cli.Configs.defaultIntValue;
import static com.hartwig.hmftools.common.cli.Configs.defaultStringValue;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.cli.Configs;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
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

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SageConfig {

    Logger LOGGER = LogManager.getLogger(SageConfig.class);

    String THREADS = "threads";
    String REFERENCE = "reference";
    String REFERENCE_BAM = "reference_bam";
    String TUMOR = "tumor";
    String TUMOR_BAM = "tumor_bam";
    String REF_GENOME = "ref_genome";
    String OUTPUT_VCF = "out";
    String MIN_MAP_QUALITY = "min_map_quality";
    String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    String PANEL_BED = "panel_bed";
    String PANEL_ONLY = "panel_only";
    String HOTSPOTS = "hotspots";
    String MAX_READ_DEPTH = "max_read_depth";
    String MAX_READ_DEPTH_PANEL = "max_read_depth_panel";
    String MAX_REALIGNMENT_DEPTH = "max_realignment_depth";
    String ASSEMBLY = "assembly";
    String CHR = "chr";
    String SLICE_SIZE = "slice_size";
    String MNV = "mnv_enabled";
    String READ_CONTEXT_FLANK_SIZE = "read_context_flank_size";

    int DEFAULT_THREADS = 2;
    int DEFAULT_MIN_MAP_QUALITY = 10;
    int DEFAULT_MAX_READ_DEPTH = 1000;
    int DEFAULT_MAX_READ_DEPTH_PANEL = 100_000;
    int DEFAULT_MAX_REALIGNMENT_DEPTH = 1000;
    int DEFAULT_SLICE_SIZE = 100_000;
    int DEFAULT_READ_CONTEXT_FLANK_SIZE = 10;
    boolean DEFAULT_MNV = true;

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(MNV, true, "Enable MNVs [" + DEFAULT_MNV + "]");
        options.addOption(ASSEMBLY, true, "Assembly, must be one of [hg19, hg38]");
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(REF_GENOME, true, "Path to indexed ref genome fasta file");
        options.addOption(OUTPUT_VCF, true, "Path to output vcf");
        options.addOption(MIN_MAP_QUALITY, true, "Min map quality to apply to non-hotspot variants [" + DEFAULT_MIN_MAP_QUALITY + "]");
        options.addOption(CHR, true, "Run for single chromosome");
        options.addOption(SLICE_SIZE, true, "Slice size [" + DEFAULT_SLICE_SIZE + "]");
        options.addOption(READ_CONTEXT_FLANK_SIZE, true, "Size of read context flank [" + DEFAULT_READ_CONTEXT_FLANK_SIZE + "]");

        options.addOption(MAX_READ_DEPTH, true, "Max depth to look for evidence [" + DEFAULT_MAX_READ_DEPTH + "]");
        options.addOption(MAX_READ_DEPTH_PANEL, true, "Max depth to look for evidence [" + DEFAULT_MAX_READ_DEPTH_PANEL + "]");
        options.addOption(MAX_REALIGNMENT_DEPTH, true, "Max depth to check for realignment [" + DEFAULT_MAX_REALIGNMENT_DEPTH + "]");
        options.addOption(HIGH_CONFIDENCE_BED, true, "High confidence regions bed file");
        options.addOption(PANEL_BED, true, "Panel regions bed file");
        options.addOption(PANEL_ONLY, false, "Only examine panel for variants");
        options.addOption(HOTSPOTS, true, "Hotspots");
        FilterConfig.createOptions().getOptions().forEach(options::addOption);
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
    List<String> tumor();

    @NotNull
    List<String> tumorBam();

    @NotNull
    String outputFile();

    @NotNull
    List<HmfTranscriptRegion> transcriptRegions();

    @NotNull
    default String baseQualityRecalibrationFile(@NotNull final String sample) {
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

    default int typicalReadLength() {
        return 151;
    }

    int regionSliceSize();

    int minMapQuality();

    int maxRealignmentDepth();

    int maxReadDepth();

    int maxReadDepthPanel();

    default int maxSkippedReferenceRegions() {
        return 50;
    }

    int readContextFlankSize();

    @NotNull
    static SageConfig createConfig(@NotNull final String version, @NotNull final CommandLine cmd) throws ParseException {
        final int threads = defaultIntValue(cmd, THREADS, DEFAULT_THREADS);
        final String assembly = cmd.getOptionValue(ASSEMBLY, "UNKNOWN");
        if (!assembly.equals("hg19") && !assembly.equals("hg38")) {
            throw new ParseException("Assembly must be one of hg19 or hg38");
        }

        final List<String> referenceList = Lists.newArrayList();
        if (cmd.hasOption(REFERENCE)) {
            referenceList.addAll(Arrays.asList(cmd.getOptionValue(REFERENCE).split(",")));
        }

        final List<String> referenceBamList = Lists.newArrayList();
        if (cmd.hasOption(REFERENCE_BAM)) {
            referenceBamList.addAll(Arrays.asList(cmd.getOptionValue(REFERENCE_BAM, Strings.EMPTY).split(",")));
        }

        if (referenceList.size() != referenceBamList.size()) {
            throw new ParseException("Each reference sample must have matching bam");
        }

        for (String referenceBam : referenceBamList) {
            if (!new File(referenceBam).exists()) {
                throw new ParseException("Unable to locate reference bam " + referenceBam);
            }
        }

        if (!cmd.hasOption(REF_GENOME)) {
            throw new ParseException(REF_GENOME + " is a mandatory argument");
        }

        if (!cmd.hasOption(OUTPUT_VCF)) {
            throw new ParseException(OUTPUT_VCF + " is a mandatory argument");
        }

        final List<String> tumorList = Lists.newArrayList();
        if (cmd.hasOption(TUMOR)) {
            tumorList.addAll(Arrays.asList(cmd.getOptionValue(TUMOR).split(",")));
        }

        final List<String> tumorBamList = Lists.newArrayList();
        if (cmd.hasOption(TUMOR_BAM)) {
            tumorBamList.addAll(Arrays.asList(cmd.getOptionValue(TUMOR_BAM, Strings.EMPTY).split(",")));
        }

        if (tumorList.size() != tumorBamList.size()) {
            throw new ParseException("Each tumor sample must have matching bam");
        }

        for (String tumorBam : tumorBamList) {
            if (!new File(tumorBam).exists()) {
                throw new ParseException("Unable to locate tumor bam " + tumorBam);
            }
        }

        if (tumorList.isEmpty()) {
            throw new ParseException("At least one tumor must be supplied");
        }

        final List<HmfTranscriptRegion> transcripts =
                assembly.equals("hg19") ? HmfGenePanelSupplier.allGeneList37() : HmfGenePanelSupplier.allGeneList38();

        final Set<String> chromosomes = Sets.newHashSet();
        final String chromosomeList = defaultStringValue(cmd, CHR, Strings.EMPTY);
        if (!chromosomeList.isEmpty()) {
            chromosomes.addAll(Lists.newArrayList(chromosomeList.split(",")));
        }

        return ImmutableSageConfig.builder()
                .version(version)
                .transcriptRegions(transcripts)
                .outputFile(cmd.getOptionValue(OUTPUT_VCF))
                .threads(threads)
                .chromosomes(chromosomes)
                .reference(referenceList)
                .referenceBam(referenceBamList)
                .tumor(tumorList)
                .tumorBam(tumorBamList)
                .mnvEnabled(defaultBooleanValue(cmd, MNV, DEFAULT_MNV))
                .refGenome(cmd.getOptionValue(REF_GENOME))
                .regionSliceSize(defaultIntValue(cmd, SLICE_SIZE, DEFAULT_SLICE_SIZE))
                .readContextFlankSize(defaultIntValue(cmd, READ_CONTEXT_FLANK_SIZE, DEFAULT_READ_CONTEXT_FLANK_SIZE))
                .minMapQuality(defaultIntValue(cmd, MIN_MAP_QUALITY, DEFAULT_MIN_MAP_QUALITY))
                .maxReadDepth(defaultIntValue(cmd, MAX_READ_DEPTH, DEFAULT_MAX_READ_DEPTH))
                .maxReadDepthPanel(defaultIntValue(cmd, MAX_READ_DEPTH_PANEL, DEFAULT_MAX_READ_DEPTH_PANEL))
                .maxRealignmentDepth(defaultIntValue(cmd, MAX_REALIGNMENT_DEPTH, DEFAULT_MAX_REALIGNMENT_DEPTH))
                .filter(FilterConfig.createConfig(cmd))
                .panelBed(cmd.getOptionValue(PANEL_BED, Strings.EMPTY))
                .highConfidenceBed(cmd.getOptionValue(HIGH_CONFIDENCE_BED, Strings.EMPTY))
                .hotspots(cmd.getOptionValue(HOTSPOTS, Strings.EMPTY))
                .qualityConfig(QualityConfig.createConfig(cmd, transcripts))
                .baseQualityRecalibrationConfig(BaseQualityRecalibrationConfig.createConfig(cmd))
                .panelOnly(Configs.containsFlag(cmd, PANEL_ONLY))
                .build();
    }
}
