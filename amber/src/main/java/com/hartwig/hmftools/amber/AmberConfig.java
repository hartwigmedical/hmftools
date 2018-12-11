package com.hartwig.hmftools.amber;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AmberConfig {

    Logger LOGGER = LogManager.getLogger(AmberConfig.class);

    int MIN_PARITION = 10000;
    int DEFAULT_THREADS = 1;
    int DEFAULT_MIN_BASE_QUALITY = 13;
    int DEFAULT_MIN_MAPPING_QUALITY = 1;
    double DEFAULT_MIN_DEPTH_PERCENTAGE = 0.5;
    double DEFAULT_MAX_DEPTH_PERCENTAGE = 1.5;
    double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    String TUMOR = "tumor";
    String BED_FILE = "bed";
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

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(BED_FILE, true, "Path to BAF locations bed file");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(MIN_BASE_QUALITY, true, "Minimum base quality for a base to be considered [" + DEFAULT_MIN_BASE_QUALITY + "]");
        options.addOption(MIN_MAPPING_QUALITY,
                true,
                "Minimum mapping quality for an alignment to be used [" + DEFAULT_MIN_MAPPING_QUALITY + "]");
        options.addOption(MIN_HET_AF_PERCENTAGE, true, "Min heterozygous AF% [" + DEFAULT_MIN_HET_AF_PERCENTAGE + "]");
        options.addOption(MAX_HET_AF_PERCENTAGE, true, "Max heterozygous AF% [" + DEFAULT_MAX_HET_AF_PERCENTAGE + "]");
        options.addOption(MIN_DEPTH_PERCENTAGE, true, "Max percentage of median depth [" + DEFAULT_MIN_DEPTH_PERCENTAGE + "]");
        options.addOption(MAX_DEPTH_PERCENTAGE, true, "Min percentage of median depth [" + DEFAULT_MAX_DEPTH_PERCENTAGE + "]");
        return options;
    }

    int threadCount();

    int minBaseQuality();

    int minMappingQuality();

    double minDepthPercent();

    double maxDepthPercent();

    double minHetAfPercent();

    double maxHetAfPercent();

    String bedFilePath();

    String tumorBamPath();

    String referenceBamPath();

    String refGenomePath();

    String outputDirectory();

    String normal();

    String tumor();

    @NotNull
    static AmberConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {
        final int threadCount = defaultIntValue(cmd, THREADS, DEFAULT_THREADS);
        final int minBaseQuality = defaultIntValue(cmd, MIN_BASE_QUALITY, DEFAULT_MIN_BASE_QUALITY);
        final int minMappingQuality = defaultIntValue(cmd, MIN_MAPPING_QUALITY, DEFAULT_MIN_MAPPING_QUALITY);

        final double minDepthPercent = defaultDoubleValue(cmd, MIN_DEPTH_PERCENTAGE, DEFAULT_MIN_DEPTH_PERCENTAGE);
        final double maxDepthPercent = defaultDoubleValue(cmd, MAX_DEPTH_PERCENTAGE, DEFAULT_MAX_DEPTH_PERCENTAGE);
        final double minHetAfPercent = defaultDoubleValue(cmd, MIN_HET_AF_PERCENTAGE, DEFAULT_MIN_HET_AF_PERCENTAGE);
        final double maxHetAfPercent = defaultDoubleValue(cmd, MAX_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);

        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String bedFilePath = parameter(cmd, BED_FILE, missingJoiner);
        final String tumorBamPath = parameter(cmd, TUMOR_BAM, missingJoiner);
        final String referenceBamPath = parameter(cmd, REFERENCE_BAM, missingJoiner);
        final String refGenomePath = parameter(cmd, REF_GENOME, missingJoiner);
        final String outputDirectory = parameter(cmd, OUTPUT_DIR, missingJoiner);
        final String normal = parameter(cmd, REFERENCE, missingJoiner);
        final String tumor = parameter(cmd, TUMOR, missingJoiner);
        final String missing = missingJoiner.toString();

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        return ImmutableAmberConfig.builder()
                .threadCount(threadCount)
                .minBaseQuality(minBaseQuality)
                .minMappingQuality(minMappingQuality)
                .minDepthPercent(minDepthPercent)
                .maxDepthPercent(maxDepthPercent)
                .minHetAfPercent(minHetAfPercent)
                .maxHetAfPercent(maxHetAfPercent)
                .bedFilePath(bedFilePath)
                .tumorBamPath(tumorBamPath)
                .referenceBamPath(referenceBamPath)
                .refGenomePath(refGenomePath)
                .outputDirectory(outputDirectory)
                .normal(normal)
                .tumor(tumor)
                .build();
    }

    static double defaultDoubleValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.valueOf(cmd.getOptionValue(opt));
            if (!Doubles.equal(result, defaultValue)) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
    }

    static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue) {
        if (cmd.hasOption(opt)) {
            final int result = Integer.valueOf(cmd.getOptionValue(opt));
            if (result != defaultValue) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
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
