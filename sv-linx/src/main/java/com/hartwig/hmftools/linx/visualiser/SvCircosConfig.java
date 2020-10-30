package com.hartwig.hmftools.linx.visualiser;

import com.hartwig.hmftools.common.utils.Doubles;

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

    String GENE_RELATIVE_SIZE = "gene_relative_size";
    String CNA_RELATIVE_SIZE = "cna_relative_size";
    String SEGMENT_RELATIVE_SIZE = "segment_relative_size";

    // Radial Arguments
    String INNER_RADIUS = "inner_radius";
    String OUTER_RADIUS = "outer_radius";
    String GAP_RADIUS = "gap_radius";
    String EXON_RANK_RADIUS = "exon_rank_radius";

    // Chromosome Range Panel Arguments
    String CHR_RANGE_HEIGHT = "chr_range_height";
    String CHR_RANGE_COLUMNS = "chr_range_columns";

    // Fusion Panel Arguments
    String FUSION_HEIGHT = "fusion_height";
    String FUSION_LEGEND_ROWS = "fusion_legend_rows";
    String FUSION_LEGEND_HEIGHT_PER_ROW = "fusion_legend_height_per_row";

    // Font Parameters
    String MIN_LABEL_SIZE = "min_label_size";
    String MAX_LABEL_SIZE = "max_label_size";
    String MAX_GENE_CHARACTERS = "max_gene_characters";
    String MAX_DISTANCE_LABELS = "max_distance_labels";
    String MAX_POSITION_LABELS = "max_position_labels";

    // Line Size
    String MIN_LINE_SIZE = "min_line_size";
    String MAX_LINE_SIZE = "max_line_size";
    String GLYPH_SIZE = "glyph_size";

    // Interpolation
    String INTERPOLATE_CNA_POSITIONS = "interpolate_cna_positions";
    String INTERPOLATE_EXON_POSITIONS = "interpolate_exon_positions";

    // Debug
    String EXACT_POSITION = "exact_position";
    String SHOW_SV_ID = "show_sv_id";
    String STEP = "step";

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
    double DEFAULT_EXON_RANK_RADIUS = 0.025;

    int DEFAULT_MIN_LABEL_SIZE = 35;
    int DEFAULT_MAX_LABEL_SIZE = 40;

    int DEFAULT_MIN_LINE_SIZE = 1;
    int DEFAULT_MAX_LINE_SIZE = 12;

    int DEFAULT_GLYPH_SIZE = 20;
    int DEFAULT_MAX_GENE_CHARACTERS = 5;

    int DEFAULT_MAX_DISTANCE_LABELS = 100;
    int DEFAULT_MAX_POSITION_LABELS = 60;

    static void addOptions(@NotNull Options options)
    {
        options.addOption(OUTER_RADIUS, true, "Outermost ending radius of chromosome track [" + DEFAULT_OUTER_RADIUS + "]");
        options.addOption(INNER_RADIUS, true, "Innermost starting radius of minor-allele ploidy track [" + DEFAULT_INNER_RADIUS + "]");
        options.addOption(GAP_RADIUS, true, "Radial gap between tracks [" + DEFAULT_GAP_RADIUS + "]");
        options.addOption(EXON_RANK_RADIUS, true, "Radial gap left for exon rank labels [" + DEFAULT_EXON_RANK_RADIUS + "]");

        options.addOption(FUSION_HEIGHT, true, "Height of each fusion in pixels [" + DEFAULT_FUSION_HEIGHT + "]");
        options.addOption(FUSION_LEGEND_ROWS, true, "Number of rows in protein domain legend  [" + DEFAULT_FUSION_LEGEND_ROWS + "]");
        options.addOption(FUSION_LEGEND_HEIGHT_PER_ROW,
                true,
                "Height of each row in protein domain legend [" + DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW + "]");

        options.addOption(CHR_RANGE_HEIGHT, true, "Chromosome range row height in pixels [" + DEFAULT_CHR_RANGE_HEIGHT + "]");
        options.addOption(CHR_RANGE_COLUMNS, true, "Chromosome range row columns [" + DEFAULT_CHR_RANGE_COLUMNS + "]");

        options.addOption(GENE_RELATIVE_SIZE,
                true,
                "Size of gene track relative to segments and copy number alterations [" + DEFAULT_GENE_RELATIVE_SIZE + "]");
        options.addOption(SEGMENT_RELATIVE_SIZE,
                true,
                "Size of segment track relative to copy number alterations and genes[" + DEFAULT_SEGMENT_RELATIVE_SIZE + "]");
        options.addOption(CNA_RELATIVE_SIZE,
                true,
                "Size of gene copy number alteration relative to genes and segments [" + DEFAULT_CNA_RELATIVE_SIZE + "]");

        options.addOption(MIN_LABEL_SIZE, true, "Minimum label size in pixels [" + DEFAULT_MIN_LABEL_SIZE + "]");
        options.addOption(MAX_LABEL_SIZE, true, "Maximum label size in pixels [" + DEFAULT_MAX_LABEL_SIZE + "]");
        options.addOption(MAX_DISTANCE_LABELS, true, "Maximum number of distance labels before removing them [" + DEFAULT_MAX_DISTANCE_LABELS + "]");
        options.addOption(MAX_POSITION_LABELS, true, "Maximum number of position labels before increasing distance between labels [" + DEFAULT_MAX_POSITION_LABELS + "]");
        options.addOption(MAX_GENE_CHARACTERS,
                true,
                "Maximum number of character in gene allowed before scaling [" + DEFAULT_MAX_GENE_CHARACTERS + "]");

        options.addOption(MIN_LINE_SIZE, true, "Minimum line size in pixels [" + DEFAULT_MIN_LINE_SIZE + "]");
        options.addOption(MAX_LINE_SIZE, true, "Maximum line size in pixels [" + DEFAULT_MAX_LINE_SIZE + "]");
        options.addOption(GLYPH_SIZE, true, "Size of glyphs in pixels [" + DEFAULT_GLYPH_SIZE + "]");

        options.addOption(INTERPOLATE_CNA_POSITIONS, false, "Interpolate copy number positions rather than adjust scale");
        options.addOption(INTERPOLATE_EXON_POSITIONS, false, "Interpolate exon positions rather than adjust scale");

        options.addOption(EXACT_POSITION, false, "Display exact positions at all break ends");
        options.addOption(SHOW_SV_ID, false, "Display SV Id next to position");
        options.addOption(STEP, false, "Create image at each step of derivative chromosome");

    }

    // ----------------------- Fusion Parameters
    int fusionLegendRows();

    int fusionLegendHeightPerRow();

    int fusionHeight();

    int chromosomeRangeHeight();

    int chromosomeRangeColumns();

    int minLabelSize();

    int maxLabelSize();

    int maxGeneCharacters();

    int maxNumberOfPositionLabels();

    int maxNumberOfDistanceLabels();

    double geneRelativeSize();

    double segmentRelativeSize();

    double copyNumberRelativeSize();

    int minLineSize();

    int maxLineSize();

    int glyphSize();

    boolean exactPosition();

    boolean showSvId();

    double outerRadius();

    double innerRadius();

    double gapRadius();

    double exonRankRadius();

    boolean interpolateCopyNumberPositions();

    boolean interpolateExonPositions();

    boolean step();

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
        if (maxLabelSize < minLabelSize)
        {
            throw new ParseException("Parameter " + MAX_LABEL_SIZE + " should be > " + MIN_LABEL_SIZE);
        }

        int minLineSize = defaultIntValue(cmd, MIN_LINE_SIZE, DEFAULT_MIN_LINE_SIZE);
        if (minLineSize <= 0)
        {
            throw new ParseException("Parameter " + MIN_LINE_SIZE + " should be > 0");
        }

        int maxLineSize = defaultIntValue(cmd, MAX_LINE_SIZE, DEFAULT_MAX_LINE_SIZE);
        if (maxLineSize < minLineSize)
        {
            throw new ParseException("Parameter " + MAX_LINE_SIZE + " should be > " + MIN_LINE_SIZE);
        }

        return ImmutableSvCircosConfig.builder()

                .outerRadius(radialParameter(cmd, OUTER_RADIUS, DEFAULT_OUTER_RADIUS))
                .innerRadius(radialParameter(cmd, INNER_RADIUS, DEFAULT_INNER_RADIUS))
                .gapRadius(radialParameter(cmd, GAP_RADIUS, DEFAULT_GAP_RADIUS))
                .exonRankRadius(radialParameter(cmd, EXON_RANK_RADIUS, DEFAULT_EXON_RANK_RADIUS))

                .geneRelativeSize(relativeParameter(cmd, GENE_RELATIVE_SIZE, DEFAULT_GENE_RELATIVE_SIZE))
                .segmentRelativeSize(relativeParameter(cmd, SEGMENT_RELATIVE_SIZE, DEFAULT_SEGMENT_RELATIVE_SIZE))
                .copyNumberRelativeSize(relativeParameter(cmd, CNA_RELATIVE_SIZE, DEFAULT_CNA_RELATIVE_SIZE))

                .minLabelSize(minLabelSize)
                .maxLabelSize(maxLabelSize)
                .maxGeneCharacters(defaultIntValue(cmd, MAX_GENE_CHARACTERS, DEFAULT_MAX_GENE_CHARACTERS))
                .maxNumberOfDistanceLabels(defaultIntValue(cmd, MAX_DISTANCE_LABELS, DEFAULT_MAX_DISTANCE_LABELS))
                .maxNumberOfPositionLabels(defaultIntValue(cmd, MAX_POSITION_LABELS, DEFAULT_MAX_POSITION_LABELS))

                .fusionLegendRows(defaultIntValue(cmd, FUSION_LEGEND_ROWS, DEFAULT_FUSION_LEGEND_ROWS))
                .fusionLegendHeightPerRow(defaultIntValue(cmd, FUSION_LEGEND_HEIGHT_PER_ROW, DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW))
                .fusionHeight(defaultIntValue(cmd, FUSION_HEIGHT, DEFAULT_FUSION_HEIGHT))

                .chromosomeRangeColumns(defaultIntValue(cmd, CHR_RANGE_COLUMNS, DEFAULT_CHR_RANGE_COLUMNS))
                .chromosomeRangeHeight(defaultIntValue(cmd, CHR_RANGE_HEIGHT, DEFAULT_CHR_RANGE_HEIGHT))

                .minLineSize(minLineSize)
                .maxLineSize(maxLineSize)
                .glyphSize(defaultIntValue(cmd, GLYPH_SIZE, DEFAULT_GLYPH_SIZE))

                .interpolateCopyNumberPositions(cmd.hasOption(INTERPOLATE_CNA_POSITIONS))
                .interpolateExonPositions(cmd.hasOption(INTERPOLATE_EXON_POSITIONS))

                .exactPosition(cmd.hasOption(EXACT_POSITION))
                .showSvId(cmd.hasOption(SHOW_SV_ID))
                .step(cmd.hasOption(STEP))

                .build();
    }

    static double relativeParameter(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue)
            throws ParseException
    {
        double value = defaultValue(cmd, opt, defaultValue);
        if (!Doubles.greaterOrEqual(value, 0))
        {
            throw new ParseException("Relative parameter " + opt + " should be > 0");
        }

        return value;
    }

    static double radialParameter(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue)
            throws ParseException
    {
        double value = defaultValue(cmd, opt, defaultValue);
        if (Doubles.lessThan(value, 0) || Doubles.greaterThan(value, 1))
        {
            throw new ParseException("Radial parameter " + opt + " should be between 0 and 1");
        }

        return value;
    }

    static double defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue)
    {
        if (cmd.hasOption(opt))
        {
            final double result = Double.parseDouble(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue)
    {
        if (cmd.hasOption(opt))
        {
            final int result = Integer.parseInt(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }
}
