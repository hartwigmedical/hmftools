package com.hartwig.hmftools.linx.visualiser;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CircosConfig
{
    public final int FusionLegendRows;
    public final int FusionLegendHeightPerRow;
    public final int FusionHeight;

    public final int ChromosomeRangeHeight;
    public final int ChromosomeRangeColumns;

    public final int MinLabelSize;
    public final int MaxLabelSize;
    public final int MaxGeneCharacters;
    public final int MaxNumberOfPositionLabels;
    public final int MaxNumberOfDistanceLabels;
    public final double GeneRelativeSize;
    public final double SegmentRelativeSize;
    public final double CopyNumberRelativeSize;

    public final int MinLineSize;
    public final int MaxLineSize;
    public final int GlyphSize;
    public final double OuterRadius;
    public final double InnerRadius;
    public final double GapRadius;
    public final double ExonRankRadius;
    public final boolean InterpolateCopyNumberPositions;
    public final boolean InterpolateExonPositions;
    public final boolean Step;

    public final boolean ExactPosition;
    public final boolean ShowSvId;

    private static final String GENE_RELATIVE_SIZE = "gene_relative_size";
    private static final String CNA_RELATIVE_SIZE = "cna_relative_size";
    private static final String SEGMENT_RELATIVE_SIZE = "segment_relative_size";

    // Radial Arguments
    private static final String INNER_RADIUS = "inner_radius";
    private static final String OUTER_RADIUS = "outer_radius";
    private static final String GAP_RADIUS = "gap_radius";
    private static final String EXON_RANK_RADIUS = "exon_rank_radius";

    // Chromosome Range Panel Arguments
    private static final String CHR_RANGE_HEIGHT = "chr_range_height";
    private static final String CHR_RANGE_COLUMNS = "chr_range_columns";

    // Fusion Panel Arguments
    private static final String FUSION_HEIGHT = "fusion_height";
    private static final String FUSION_LEGEND_ROWS = "fusion_legend_rows";
    private static final String FUSION_LEGEND_HEIGHT_PER_ROW = "fusion_legend_height_per_row";

    // Font Parameters
    private static final String MIN_LABEL_SIZE = "min_label_size";
    private static final String MAX_LABEL_SIZE = "max_label_size";
    private static final String MAX_GENE_CHARACTERS = "max_gene_characters";
    private static final String MAX_DISTANCE_LABELS = "max_distance_labels";
    private static final String MAX_POSITION_LABELS = "max_position_labels";

    // Line Size
    private static final String MIN_LINE_SIZE = "min_line_size";
    private static final String MAX_LINE_SIZE = "max_line_size";
    private static final String GLYPH_SIZE = "glyph_size";

    // Interpolation
    private static final String INTERPOLATE_CNA_POSITIONS = "interpolate_cna_positions";
    private static final String INTERPOLATE_EXON_POSITIONS = "interpolate_exon_positions";

    // Debug
    private static final String EXACT_POSITION = "exact_position";
    private static final String SHOW_SV_ID = "show_sv_id";
    private static final String STEP = "step";

    private static final int DEFAULT_FUSION_HEIGHT = 250;
    private static final int DEFAULT_FUSION_LEGEND_ROWS = 1;
    private static final int DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW = 35;

    private static final int DEFAULT_CHR_RANGE_HEIGHT = 150;
    private static final int DEFAULT_CHR_RANGE_COLUMNS = 6;

    private static final double DEFAULT_GENE_RELATIVE_SIZE = 0.3;
    private static final double DEFAULT_SEGMENT_RELATIVE_SIZE = 1;
    private static final double DEFAULT_CNA_RELATIVE_SIZE = 2;

    private static final double DEFAULT_OUTER_RADIUS = 0.88;
    private static final double DEFAULT_INNER_RADIUS = 0.20;
    private static final double DEFAULT_GAP_RADIUS = 0.025;
    private static final double DEFAULT_EXON_RANK_RADIUS = 0.025;

    private static final int DEFAULT_MIN_LABEL_SIZE = 35;
    private static final int DEFAULT_MAX_LABEL_SIZE = 40;

    private static final int DEFAULT_MIN_LINE_SIZE = 1;
    private static final int DEFAULT_MAX_LINE_SIZE = 12;

    private static final int DEFAULT_GLYPH_SIZE = 20;
    private static final int DEFAULT_MAX_GENE_CHARACTERS = 5;

    private static final int DEFAULT_MAX_DISTANCE_LABELS = 100;
    private static final int DEFAULT_MAX_POSITION_LABELS = 60;

    public static void addOptions(@NotNull Options options)
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

    public long labelSize(long count)
    {
        if (count > MaxNumberOfDistanceLabels)
        {
            return MinLabelSize;
        }

        return Math.round(MaxLabelSize - 1d * count * (MaxLabelSize - MinLabelSize) / MaxNumberOfDistanceLabels);
    }

    public CircosConfig(@NotNull final CommandLine cmd) throws ParseException
    {
        MinLabelSize = getConfigValue(cmd, MIN_LABEL_SIZE, DEFAULT_MIN_LABEL_SIZE);
        MaxLabelSize = getConfigValue(cmd, MAX_LABEL_SIZE, DEFAULT_MAX_LABEL_SIZE);

        OuterRadius = getConfigValue(cmd, OUTER_RADIUS, DEFAULT_OUTER_RADIUS);
        InnerRadius = getConfigValue(cmd, INNER_RADIUS, DEFAULT_INNER_RADIUS);
        GapRadius = getConfigValue(cmd, GAP_RADIUS, DEFAULT_GAP_RADIUS);
        ExonRankRadius = getConfigValue(cmd, EXON_RANK_RADIUS, DEFAULT_EXON_RANK_RADIUS);

        GeneRelativeSize = getConfigValue(cmd, GENE_RELATIVE_SIZE, DEFAULT_GENE_RELATIVE_SIZE);
        SegmentRelativeSize = getConfigValue(cmd, SEGMENT_RELATIVE_SIZE, DEFAULT_SEGMENT_RELATIVE_SIZE);
        CopyNumberRelativeSize = getConfigValue(cmd, CNA_RELATIVE_SIZE, DEFAULT_CNA_RELATIVE_SIZE);

        MaxGeneCharacters = getConfigValue(cmd, MAX_GENE_CHARACTERS, DEFAULT_MAX_GENE_CHARACTERS);
        MaxNumberOfDistanceLabels = getConfigValue(cmd, MAX_DISTANCE_LABELS, DEFAULT_MAX_DISTANCE_LABELS);
        MaxNumberOfPositionLabels = getConfigValue(cmd, MAX_POSITION_LABELS, DEFAULT_MAX_POSITION_LABELS);

        FusionLegendRows = getConfigValue(cmd, FUSION_LEGEND_ROWS, DEFAULT_FUSION_LEGEND_ROWS);
        FusionLegendHeightPerRow = getConfigValue(cmd, FUSION_LEGEND_HEIGHT_PER_ROW, DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW);
        FusionHeight = getConfigValue(cmd, FUSION_HEIGHT, DEFAULT_FUSION_HEIGHT);

        ChromosomeRangeColumns = getConfigValue(cmd, CHR_RANGE_COLUMNS, DEFAULT_CHR_RANGE_COLUMNS);
        ChromosomeRangeHeight = getConfigValue(cmd, CHR_RANGE_HEIGHT, DEFAULT_CHR_RANGE_HEIGHT);

        MinLineSize = getConfigValue(cmd, MIN_LINE_SIZE, DEFAULT_MIN_LINE_SIZE);
        MaxLineSize = getConfigValue(cmd, MAX_LINE_SIZE, DEFAULT_MAX_LINE_SIZE);
        GlyphSize = getConfigValue(cmd, GLYPH_SIZE, DEFAULT_GLYPH_SIZE);

        InterpolateCopyNumberPositions = cmd.hasOption(INTERPOLATE_CNA_POSITIONS);
        InterpolateExonPositions = cmd.hasOption(INTERPOLATE_EXON_POSITIONS);

        ExactPosition = cmd.hasOption(EXACT_POSITION);
        ShowSvId = cmd.hasOption(SHOW_SV_ID);
        Step = cmd.hasOption(STEP);
    }

    public boolean isValid()
    {
        if (MinLabelSize <= 0)
        {
            VIS_LOGGER.error("Parameter " + MIN_LABEL_SIZE + " should be > 0");
            return false;
        }

        if (MaxLabelSize < MinLabelSize)
        {
            VIS_LOGGER.error("Parameter " + MAX_LABEL_SIZE + " should be > " + MIN_LABEL_SIZE);
            return false;
        }

        if (MinLineSize <= 0)
        {
            VIS_LOGGER.error("Parameter " + MIN_LINE_SIZE + " should be > 0");
            return false;
        }

        if (MaxLineSize < MinLineSize)
        {
            VIS_LOGGER.error("Parameter " + MAX_LINE_SIZE + " should be > " + MIN_LINE_SIZE);
            return false;
        }

        if(!validRadialParameter(OUTER_RADIUS, OuterRadius))
            return false;

        if(!validRadialParameter(INNER_RADIUS, InnerRadius))
            return false;

        if(!validRadialParameter(GAP_RADIUS, GapRadius))
            return false;

        if(!validRadialParameter(EXON_RANK_RADIUS, ExonRankRadius))
            return false;

        if(!validRelativeParameter(GENE_RELATIVE_SIZE, GeneRelativeSize))
            return false;

        if(!validRelativeParameter(SEGMENT_RELATIVE_SIZE, SegmentRelativeSize))
            return false;

        if(!validRelativeParameter(CNA_RELATIVE_SIZE, CopyNumberRelativeSize))
            return false;

        return true;
    }

    private static boolean validRelativeParameter(final String opt, double value)
    {
        if(!Doubles.greaterOrEqual(value, 0))
        {
            VIS_LOGGER.error("relative parameter({}) value({}), must be greater than 0", opt, value);
            return false;
        }

        return true;
    }

    private static boolean validRadialParameter(final String opt, double value)
    {
        if(Doubles.lessThan(value, 0) || Doubles.greaterThan(value, 1))
        {
            VIS_LOGGER.error("radial parameter({}) value({}), must be betwee 0 - 1", opt, value);
            return false;
        }

        return true;
    }
}
