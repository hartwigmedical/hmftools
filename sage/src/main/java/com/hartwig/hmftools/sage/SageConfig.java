package com.hartwig.hmftools.sage;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

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
public interface SageConfig {

    Logger LOGGER = LogManager.getLogger(SageConfig.class);

    String THREADS = "threads";
    String REFERENCE = "reference";
    String REFERENCE_BAM = "reference_bam";
    String TUMOR = "tumor";
    String TUMOR_BAM = "tumor_bam";
    String REF_GENOME = "ref_genome";
    String OUTPUT_VCF = "out";

    int DEFAULT_THREADS = 2;

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(REF_GENOME, true, "Path to indexed ref genome fasta file");
        options.addOption(OUTPUT_VCF, true, "Path to output vcf");

        return options;
    }

    int threads();

    @NotNull
    String reference();

    @NotNull
    String referenceBam();

    @NotNull
    String refGenome();

    @NotNull
    List<String> tumor();

    @NotNull
    List<String> tumorBam();

    @NotNull
    String outputFile();

    default int minQuality() {
        return 13;
    }

    @NotNull
    static SageConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {

        final int threads = defaultIntValue(cmd, THREADS, DEFAULT_THREADS);
        final String reference = cmd.getOptionValue(REFERENCE);
        final String reference_bam = cmd.getOptionValue(REFERENCE_BAM);

        final List<String> tumorList = Lists.newArrayList();
        tumorList.addAll(Arrays.asList(cmd.getOptionValue(TUMOR).split(",")));

        final List<String> tumorBamList = Lists.newArrayList();
        tumorBamList.addAll(Arrays.asList(cmd.getOptionValue(TUMOR_BAM).split(",")));

        if (tumorList.size() != tumorBamList.size()) {
            throw new ParseException("TODO");
        }

        return ImmutableSageConfig.builder()
                .outputFile(cmd.getOptionValue(OUTPUT_VCF))
                .threads(threads)
                .reference(reference)
                .referenceBam(reference_bam)
                .tumor(tumorList)
                .tumorBam(tumorBamList)
                .refGenome(cmd.getOptionValue(REF_GENOME))
                .build();

    }

    @NotNull
    static String defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, @NotNull final String defaultValue) {
        return cmd.hasOption(opt) ? cmd.getOptionValue(opt) : defaultValue;
    }

    static double defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue) {
        if (cmd.hasOption(opt)) {
            final int result = Integer.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

}
