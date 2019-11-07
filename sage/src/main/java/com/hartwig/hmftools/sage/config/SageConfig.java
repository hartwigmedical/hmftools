package com.hartwig.hmftools.sage.config;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

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
    String MIN_BASE_QUALITY = "min_base_quality";
    String PANEL = "panel";
    String HOTSPOTS = "hotspots";

    int DEFAULT_THREADS = 2;
    int DEFAULT_MIN_MAP_QUALITY = 0;
    int DEFAULT_MIN_BASE_QUALITY = 13;

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
        options.addOption(MIN_MAP_QUALITY, true, "Min map quality [" + DEFAULT_MIN_MAP_QUALITY + "]");
        options.addOption(MIN_BASE_QUALITY, true, "Min base quality [" + DEFAULT_MIN_BASE_QUALITY + "]");

        options.addOption(PANEL, true, "Panel");
        options.addOption(HOTSPOTS, true, "Hotspots");
        FilterConfig.createOptions().getOptions().forEach(options::addOption);
        QualityConfig.createOptions().getOptions().forEach(options::addOption);

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

    @NotNull
    String panel();

    @NotNull
    String hotspots();

    @NotNull
    FilterConfig filter();

    @NotNull
    QualityConfig qualityConfig();

    default int regionSliceSize() {
        return 1_000_000;
    }

    int minMapQuality();

    int minBaseQuality();

    default int maxDepthCoverage() {
        return 1000;
    }

    default int maxNormalAltSupport() {
        return 3;
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
                .minMapQuality(defaultIntValue(cmd, MIN_MAP_QUALITY, DEFAULT_MIN_MAP_QUALITY))
                .minBaseQuality(defaultIntValue(cmd, MIN_BASE_QUALITY, DEFAULT_MIN_BASE_QUALITY))
                .filter(FilterConfig.createConfig(cmd))
                .panel(cmd.getOptionValue(PANEL, Strings.EMPTY))
                .hotspots(cmd.getOptionValue(HOTSPOTS, Strings.EMPTY))
                .qualityConfig(QualityConfig.createConfig(cmd))
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
