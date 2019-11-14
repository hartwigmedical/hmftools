package com.hartwig.hmftools.cobalt;

import java.util.StringJoiner;

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
public interface CobaltConfig {
    Logger LOGGER = LogManager.getLogger(CobaltConfig.class);

    int DEFAULT_THREADS = 4;
    int DEFAULT_MIN_MAPPING_QUALITY = 10;

    String THREADS = "threads";
    String REFERENCE = "reference";
    String REFERENCE_BAM = "reference_bam";
    String TUMOR = "tumor";
    String TUMOR_BAM = "tumor_bam";
    String REF_GENOME = "ref_genome";
    String OUTPUT_DIR = "output_dir";
    String GC_PROFILE = "gc_profile";
    String MIN_MAPPING_QUALITY = "min_quality";

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(MIN_MAPPING_QUALITY, true, "Min quality [" + DEFAULT_MIN_MAPPING_QUALITY + "]");
        options.addOption(GC_PROFILE, true, "Location of GC Profile");
        options.addOption(REF_GENOME, true, "Path to reference genome fasta file if using CRAM files");

        return options;
    }

    int threadCount();

    int minMappingQuality();

    @NotNull
    String gcProfilePath();

    @NotNull
    String tumorBamPath();

    @NotNull
    String referenceBamPath();

    @NotNull
    String refGenomePath();

    @NotNull
    String outputDirectory();

    @NotNull
    String reference();

    @NotNull
    String tumor();

    default int windowSize() {
        return 1000;
    }

    @NotNull
    static CobaltConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {
        final int threadCount = defaultIntValue(cmd, THREADS, DEFAULT_THREADS);
        final int minMappingQuality = defaultIntValue(cmd, MIN_MAPPING_QUALITY, DEFAULT_MIN_MAPPING_QUALITY);
        final String refGenomePath = cmd.getOptionValue(REF_GENOME, "");

        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String gcProfilePath = parameter(cmd, GC_PROFILE, missingJoiner);
        final String tumorBamPath = parameter(cmd, TUMOR_BAM, missingJoiner);
        final String referenceBamPath = parameter(cmd, REFERENCE_BAM, missingJoiner);
        final String outputDirectory = parameter(cmd, OUTPUT_DIR, missingJoiner);
        final String normal = parameter(cmd, REFERENCE, missingJoiner);
        final String tumor = parameter(cmd, TUMOR, missingJoiner);
        final String missing = missingJoiner.toString();

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        return ImmutableCobaltConfig.builder()
                .threadCount(threadCount)
                .minMappingQuality(minMappingQuality)
                .gcProfilePath(gcProfilePath)
                .tumorBamPath(tumorBamPath)
                .referenceBamPath(referenceBamPath)
                .refGenomePath(refGenomePath)
                .outputDirectory(outputDirectory)
                .reference(normal)
                .tumor(tumor)
                .build();
    }

    static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue) {
        if (cmd.hasOption(opt)) {
            final int result = Integer.parseInt(cmd.getOptionValue(opt));
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
