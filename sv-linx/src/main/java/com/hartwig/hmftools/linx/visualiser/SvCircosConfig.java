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
    String SEGMENT_RELATIVE_SIZE = "segment_relative_size";

    String MIN_LABEL_SIZE = "min_label_size";
    String MAX_LABEL_SIZE = "max_label_size";
    String OUTER_RADIUS = "outer_radius";
    String INNER_RADIUS = "inner_radius";
    String EXACT_POSITION = "exact_position";
    String SHOW_SV_ID = "show_sv_id";

    String CHR_RANGE_HEIGHT = "chr_range_height";
    String CHR_RANGE_COLUMNS = "chr_range_columns";

    String FUSION_HEIGHT = "fusion_height";
    String FUSION_LEGEND_ROWS = "fusion_legend_rows";
    String FUSION_LEGEND_HEIGHT_PER_ROW = "fusion_legend_height_per_row";

    String INTERPOLATE_CNA_POSITIONS = "interpolate_cna_positions";

    int DEFAULT_FUSION_HEIGHT = 250;
    int DEFAULT_FUSION_LEGEND_ROWS = 1;
    int DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW = 35;

    int DEFAULT_CHR_RANGE_HEIGHT = 150;
    int DEFAULT_CHR_RANGE_COLUMNS = 6;

    double DEFAULT_GENE_RELATIVE_SIZE = 0.3;
    double DEFAULT_SEGMENT_RELATIVE_SIZE = 1;
    double DEFAULT_CNA_RELATIVE_SIZE = 2;

    double DEFAULT_OUTER_RADIUS = 0.88;
    double DEFAULT_INNER_RADIUS = 0.20;
    double DEFAULT_GAP_RADIUS = 0.025;

    int DEFAULT_MIN_LABEL_SIZE = 35;
    int DEFAULT_MAX_LABEL_SIZE = 40;

    static void addOptions(@NotNull Options options)
    {

        options.addOption(FUSION_HEIGHT, true, "Height of each fusion in pixels [" + DEFAULT_FUSION_HEIGHT + "]");
        options.addOption(FUSION_LEGEND_ROWS, true, "Number of rows in protein domain legend  [" + DEFAULT_FUSION_LEGEND_ROWS + "]");
        options.addOption(FUSION_LEGEND_HEIGHT_PER_ROW, true,
                "Height of each row in protein domain legend [" + DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW
                        + "]");

        options.addOption(CHR_RANGE_HEIGHT, true, "Chromosome range row height in pixels [" + DEFAULT_CHR_RANGE_HEIGHT + "]");
        options.addOption(CHR_RANGE_COLUMNS, true, "Chromosome range row columns [" + DEFAULT_CHR_RANGE_COLUMNS + "]");

        options.addOption(EXACT_POSITION, false, "Display exact position of structural variants");
        options.addOption(SHOW_SV_ID, false, "Display SV Id next to position");

        options.addOption(GENE_RELATIVE_SIZE, true,
                "Size of gene track relative to segments and copy number alterations [" + DEFAULT_GENE_RELATIVE_SIZE + "]");
        options.addOption(SEGMENT_RELATIVE_SIZE, true,
                "Size of segment track relative to copy number alterations and genes[" + DEFAULT_SEGMENT_RELATIVE_SIZE + "]");
        options.addOption(CNA_RELATIVE_SIZE, true,
                "Size of gene copy number alteration relative to genes and segments [" + DEFAULT_CNA_RELATIVE_SIZE + "]");

        options.addOption(OUTER_RADIUS, true, "Set outer radius [" + DEFAULT_OUTER_RADIUS + "]");
        options.addOption(INNER_RADIUS, true, "Set inner radius [" + DEFAULT_INNER_RADIUS + "]");

        options.addOption(MIN_LABEL_SIZE, true, "Minimum label size [" + DEFAULT_MIN_LABEL_SIZE + "]");
        options.addOption(MAX_LABEL_SIZE, true, "Maximum label size [" + DEFAULT_MAX_LABEL_SIZE + "]");
        options.addOption(INTERPOLATE_CNA_POSITIONS, false, "Interpolate copy number positions rather than adjust scale");
    }

    // ----------------------- Fusion Parameters
    int fusionLegendRows();

    int fusionLegendHeightPerRow();

    int fusionHeight();

    // ----------------------- Chromosome Parameters
    int chromosomeRangeHeight();

    int chromosomeRangeColumns();

    // ----------------------- Label Size Parameters
    int minLabelSize();

    int maxLabelSize();

    int maxLabelCharacters();

    default int maxNumberOfDistanceLabels()
    {
        return 100;
    }

    // ----------------------- Relative Size Parameters
    double geneRelativeSize();

    double segmentRelativeSize();

    double copyNumberRelativeSize();

    // ----------------------- Other
    boolean exactPosition();

    boolean showSvId();

    double outerRadius();

    double innerRadius();

    boolean interpolateCopyNumberPositions();

    default double gapRadius()
    {
        return DEFAULT_GAP_RADIUS;
    }

    default long labelSize(long count)
    {
        if (count > maxNumberOfDistanceLabels())
        {
            return minLabelSize();
        }

        return Math.round(maxLabelSize() - 1d * count * (maxLabelSize() - minLabelSize()) / maxNumberOfDistanceLabels());
    }

    @NotNull
    static SvCircosConfig createConfig(@NotNull final CommandLine cmd) throws ParseException
    {
        int minLabelSize = defaultIntValue(cmd, MIN_LABEL_SIZE, DEFAULT_MIN_LABEL_SIZE);
        if (minLabelSize <= 0)
        {
            throw new ParseException("Parameter " + MIN_LABEL_SIZE + " should be > 0");
        }

        int maxLabelSize = defaultIntValue(cmd, MAX_LABEL_SIZE, DEFAULT_MAX_LABEL_SIZE);
        if (maxLabelSize <= 0)
        {
            throw new ParseException("Parameter " + MAX_LABEL_SIZE + " should be > 0");
        }

        // TODO: Add validation on ideogram radius etc

        return ImmutableSvCircosConfig.builder()
                .fusionLegendRows(defaultIntValue(cmd, FUSION_LEGEND_ROWS, DEFAULT_FUSION_LEGEND_ROWS))
                .fusionLegendHeightPerRow(defaultIntValue(cmd, FUSION_LEGEND_HEIGHT_PER_ROW, DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW))
                .fusionHeight(defaultIntValue(cmd, FUSION_HEIGHT, DEFAULT_FUSION_HEIGHT))
                .chromosomeRangeColumns(defaultIntValue(cmd, CHR_RANGE_COLUMNS, DEFAULT_CHR_RANGE_COLUMNS))
                .chromosomeRangeHeight(defaultIntValue(cmd, CHR_RANGE_HEIGHT, DEFAULT_CHR_RANGE_HEIGHT))
                .exactPosition(cmd.hasOption(EXACT_POSITION))
                .showSvId(cmd.hasOption(SHOW_SV_ID))
                .outerRadius(defaultValue(cmd, OUTER_RADIUS, DEFAULT_OUTER_RADIUS))
                .innerRadius(defaultValue(cmd, INNER_RADIUS, DEFAULT_INNER_RADIUS))
                .geneRelativeSize(defaultValue(cmd, GENE_RELATIVE_SIZE, DEFAULT_GENE_RELATIVE_SIZE))
                .segmentRelativeSize(defaultValue(cmd, SEGMENT_RELATIVE_SIZE, DEFAULT_SEGMENT_RELATIVE_SIZE))
                .copyNumberRelativeSize(defaultValue(cmd, CNA_RELATIVE_SIZE, DEFAULT_CNA_RELATIVE_SIZE))
                .interpolateCopyNumberPositions(cmd.hasOption(INTERPOLATE_CNA_POSITIONS))
                .minLabelSize(minLabelSize)
                .maxLabelSize(maxLabelSize)
                .maxLabelCharacters(5)
                .build();
    }

    static double defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue)
    {
        if (cmd.hasOption(opt))
        {
            final double result = Double.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue)
    {
        if (cmd.hasOption(opt))
        {
            final int result = Integer.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

}
