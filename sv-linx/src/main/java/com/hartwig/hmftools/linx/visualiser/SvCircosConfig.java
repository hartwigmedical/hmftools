package com.hartwig.hmftools.linx.visualiser;

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
public interface SvCircosConfig
{
    Logger LOGGER = LogManager.getLogger(SvCircosConfig.class);

    String DISPLAY_POSITION = "displayPosition";
    String GENE_RELATIVE_SIZE = "gene_relative_size";
    String CNA_RELATIVE_SIZE = "cna_relative_size";

    String MIN_LABEL_SIZE = "min_label_size";
    String MAX_LABEL_SIZE = "max_label_size";
    String IDEOGRAM_RADIUS = "ideogram_radius";
    String EXACT_POSITION = "exact_position";

    double DEFAULT_GENE_RELATIVE_SIZE = 0.3;
    double DEFAULT_SEGMENT_RELATIVE_SIZE = 1;
    double DEFAULT_CNA_RELATIVE_SIZE = 2;

    double DEFAULT_IDEOGRAM_RADIUS = 0.88;
    double DEFAULT_INNER_RADIUS = 0.2;
    double DEFAULT_GAP_RADIUS = 0.025;

    int DEFAULT_MIN_LABEL_SIZE = 35;
    int DEFAULT_MAX_LABEL_SIZE = 40;

    boolean exactPosition();

    double geneRelativeSize();

    double segmentRelativeSize();

    double copyNumberRelativeSize();

    double ideogramRadius();

    default double gapRadius()
    {
        return DEFAULT_GAP_RADIUS;
    }

    default double innerRadius()
    {
        return DEFAULT_INNER_RADIUS;
    }

    int minLabelSize();

    int maxLabelSize();

    default int maxDistanceLabels()
    {
        return 100;
    }

    default long labelSize(long count)
    {
        if (count > maxDistanceLabels())
        {
            return minLabelSize();
        }

        return Math.round(maxLabelSize() - 1d * count * (maxLabelSize() - minLabelSize()) / maxDistanceLabels());
    }

    static void addOptions(@NotNull Options options)
    {
        options.addOption(EXACT_POSITION, false, "Display exact position of structural variants");
        options.addOption(IDEOGRAM_RADIUS, true, "Set ideogram radius [" + DEFAULT_IDEOGRAM_RADIUS + "]");
        options.addOption(MIN_LABEL_SIZE, true, "Set ideogram radius [" + DEFAULT_MIN_LABEL_SIZE + "]");
        options.addOption(MAX_LABEL_SIZE, true, "Set ideogram radius [" + DEFAULT_MAX_LABEL_SIZE + "]");
    }

    @NotNull
    static SvCircosConfig createConfig(@NotNull final CommandLine cmd) throws ParseException
    {
        int minLabelSize = defaultIntValue(cmd, MIN_LABEL_SIZE, DEFAULT_MIN_LABEL_SIZE);
        if (minLabelSize <= 0) {
            throw new ParseException("Parameter " + MIN_LABEL_SIZE + " should be > 0");
        }

        int maxLabelSize = defaultIntValue(cmd, MAX_LABEL_SIZE, DEFAULT_MAX_LABEL_SIZE);
        if (maxLabelSize <= 0) {
            throw new ParseException("Parameter " + MAX_LABEL_SIZE + " should be > 0");
        }

        // TODO: Add validation on ideogram radius etc

        return ImmutableSvCircosConfig.builder()
                .exactPosition(cmd.hasOption(EXACT_POSITION))
                .ideogramRadius(defaultValue(cmd, IDEOGRAM_RADIUS, DEFAULT_IDEOGRAM_RADIUS))
                .geneRelativeSize(DEFAULT_GENE_RELATIVE_SIZE)
                .segmentRelativeSize(DEFAULT_SEGMENT_RELATIVE_SIZE)
                .copyNumberRelativeSize(DEFAULT_CNA_RELATIVE_SIZE)
                .minLabelSize(minLabelSize)
                .maxLabelSize(maxLabelSize)
                .build();
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
