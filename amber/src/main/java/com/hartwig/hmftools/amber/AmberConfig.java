package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.utils.ConfigUtils.defaultEnumValue;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.ValidationStringency;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AmberConfig {

    int DEFAULT_THREADS = 1;
    int DEFAULT_MIN_BASE_QUALITY = 13;
    int DEFAULT_MIN_PARTITION = 10000;
    int DEFAULT_MIN_MAPPING_QUALITY = 1;
    int DEFAULT_TYPICAL_READ_DEPTH = 151;
    int DEFAULT_TUMOR_ONLY_MIN_SUPPORT = 2;
    double DEFAULT_TUMOR_ONLY_MIN_VAF = 0.05;
    double DEFAULT_MIN_DEPTH_PERCENTAGE = 0.5;
    double DEFAULT_MAX_DEPTH_PERCENTAGE = 1.5;
    double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    String TUMOR = "tumor";
    String BAF_LOCI = "loci";
    String THREADS = "threads";
    String REFERENCE = "reference";
    String TUMOR_BAM = "tumor_bam";
    String REF_GENOME = "ref_genome";
    String OUTPUT_DIR = "output_dir";
    String REFERENCE_BAM = "reference_bam";
    String MIN_BASE_QUALITY = "min_base_quality";
    String MIN_MAPPING_QUALITY = "min_mapping_quality";
    String MIN_DEPTH_PERCENTAGE = "min_depth_percent";
    String MAX_DEPTH_PERCENTAGE = "max_depth_percent";
    String MIN_HET_AF_PERCENTAGE = "min_het_af_percent";
    String MAX_HET_AF_PERCENTAGE = "max_het_af_percent";
    String VALIDATION_STRINGENCY = "validation_stringency";

    String TUMOR_ONLY = "tumor_only";
    String TUMOR_ONLY_MIN_VAF = "tumor_only_min_vaf";
    String TUMOR_ONLY_MIN_SUPPORT = "tumor_only_min_support";

    @NotNull
    static org.apache.commons.cli.Options createOptions() {
        final org.apache.commons.cli.Options options = new org.apache.commons.cli.Options();
        options.addOption(TUMOR_ONLY, false, "Tumor only mode");
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(BAF_LOCI, true, "Path to BAF loci vcf file");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(MIN_BASE_QUALITY, true, "Minimum quality for a base to be considered [" + DEFAULT_MIN_BASE_QUALITY + "]");
        options.addOption(MIN_MAPPING_QUALITY,
                true,
                "Minimum mapping quality for an alignment to be used [" + DEFAULT_MIN_MAPPING_QUALITY + "]");
        options.addOption(MIN_HET_AF_PERCENTAGE, true, "Min heterozygous AF% [" + DEFAULT_MIN_HET_AF_PERCENTAGE + "]");
        options.addOption(MAX_HET_AF_PERCENTAGE, true, "Max heterozygous AF% [" + DEFAULT_MAX_HET_AF_PERCENTAGE + "]");
        options.addOption(MIN_DEPTH_PERCENTAGE, true, "Max percentage of median depth [" + DEFAULT_MIN_DEPTH_PERCENTAGE + "]");
        options.addOption(MAX_DEPTH_PERCENTAGE, true, "Min percentage of median depth [" + DEFAULT_MAX_DEPTH_PERCENTAGE + "]");

        options.addOption(TUMOR_ONLY_MIN_VAF, true, "Min support in ref and alt in tumor only mode [" + DEFAULT_TUMOR_ONLY_MIN_VAF + "]");
        options.addOption(TUMOR_ONLY_MIN_SUPPORT,
                true,
                "Min VAF in ref and alt in tumor only mode [" + DEFAULT_TUMOR_ONLY_MIN_SUPPORT + "]");
        options.addOption(VALIDATION_STRINGENCY, true, "SAM validation strategy: STRICT, SILENT, LENIENT [STRICT]");

        return options;
    }

    boolean tumorOnly();

    int tumorOnlyMinSupport();

    double tumorOnlyMinVaf();

    int threadCount();

    int minBaseQuality();

    int minMappingQuality();

    double minDepthPercent();

    double maxDepthPercent();

    double minHetAfPercent();

    double maxHetAfPercent();

    @NotNull
    String bafLociPath();

    @NotNull
    String tumorBamPath();

    @NotNull
    List<String> reference();

    @NotNull
    List<String> referenceBamPath();

    @NotNull
    String refGenomePath();

    @NotNull
    String outputDirectory();

    @NotNull
    String tumor();

    @NotNull
    default String primaryReference() {
        return reference().get(0);
    }

    @NotNull
    default List<String> allSamples() {
        List<String> samples = Lists.newArrayList();
        samples.addAll(reference());
        samples.add(tumor());
        return samples;
    }

    @NotNull
    ValidationStringency validationStringency();

    default int typicalReadDepth() {
        return DEFAULT_TYPICAL_READ_DEPTH;
    }

    default int minPartition() {
        return DEFAULT_MIN_PARTITION;
    }

    @NotNull
    static AmberConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {
        final boolean isTumorOnly = cmd.hasOption(TUMOR_ONLY);

        final int threadCount = getConfigValue(cmd, THREADS, DEFAULT_THREADS);
        final int minBaseQuality = getConfigValue(cmd, MIN_BASE_QUALITY, DEFAULT_MIN_BASE_QUALITY);
        final int minMappingQuality = getConfigValue(cmd, MIN_MAPPING_QUALITY, DEFAULT_MIN_MAPPING_QUALITY);
        final int tumorOnlyMinSupport = getConfigValue(cmd, TUMOR_ONLY_MIN_SUPPORT, DEFAULT_TUMOR_ONLY_MIN_SUPPORT);

        final double minDepthPercent = getConfigValue(cmd, MIN_DEPTH_PERCENTAGE, DEFAULT_MIN_DEPTH_PERCENTAGE);
        final double maxDepthPercent = getConfigValue(cmd, MAX_DEPTH_PERCENTAGE, DEFAULT_MAX_DEPTH_PERCENTAGE);
        final double minHetAfPercent = getConfigValue(cmd, MIN_HET_AF_PERCENTAGE, DEFAULT_MIN_HET_AF_PERCENTAGE);
        final double maxHetAfPercent = getConfigValue(cmd, MAX_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);
        final double tumorOnlyMinVaf = getConfigValue(cmd, TUMOR_ONLY_MIN_VAF, DEFAULT_TUMOR_ONLY_MIN_VAF);

        final StringJoiner missingJoiner = new StringJoiner(", ");
        final List<String> reference =
                isTumorOnly ? Collections.emptyList() : Arrays.asList(parameter(cmd, REFERENCE, missingJoiner).split(","));
        final List<String> referenceBamPath =
                isTumorOnly ? Collections.emptyList() : Arrays.asList(parameter(cmd, REFERENCE_BAM, missingJoiner).split(","));
        final String bafLociPath = parameter(cmd, BAF_LOCI, missingJoiner);
        final String tumorBamPath = parameter(cmd, TUMOR_BAM, missingJoiner);
        final String outputDirectory = parameter(cmd, OUTPUT_DIR, missingJoiner);
        final String tumor = parameter(cmd, TUMOR, missingJoiner);
        final String missing = missingJoiner.toString();
        final String refGenomePath = cmd.getOptionValue(REF_GENOME, Strings.EMPTY);

        if (reference.size() != referenceBamPath.size()) {
            throw new ParseException("Each reference sample must have matching bam");
        }

        final ValidationStringency validationStringency =
                defaultEnumValue(cmd, VALIDATION_STRINGENCY, ValidationStringency.DEFAULT_STRINGENCY);

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        return ImmutableAmberConfig.builder()
                .tumorOnly(isTumorOnly)
                .tumorOnlyMinVaf(tumorOnlyMinVaf)
                .tumorOnlyMinSupport(tumorOnlyMinSupport)
                .threadCount(threadCount)
                .minBaseQuality(minBaseQuality)
                .minMappingQuality(minMappingQuality)
                .minDepthPercent(minDepthPercent)
                .maxDepthPercent(maxDepthPercent)
                .minHetAfPercent(minHetAfPercent)
                .maxHetAfPercent(maxHetAfPercent)
                .bafLociPath(bafLociPath)
                .tumorBamPath(tumorBamPath)
                .referenceBamPath(referenceBamPath)
                .refGenomePath(refGenomePath)
                .outputDirectory(outputDirectory)
                .reference(reference)
                .tumor(tumor)
                .validationStringency(validationStringency)
                .build();
    }

    @NotNull
    static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing) {
        final String value = cmd.getOptionValue(parameter);
        if (value == null) {
            missing.add(parameter);
            return "";
        }
        return value;
    }
}
