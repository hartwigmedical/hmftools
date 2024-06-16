package com.hartwig.hmftools.linx.visualiser;

import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class CircosConfig
{
    public final int FusionLegendRows;
    public final int FusionLegendHeightPerRow;
    public final int FusionHeight;

    public final int ChromosomeRangeHeight;
    public final int ChromosomeRangeColumns;

    public final int MaxPlotSvCount;

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

    private static final String MAX_PLOT_SVS = "max_plot_svs";

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

    private static final int DEFAULT_MAX_PLOT_SV_COUNT = 2000;

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

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(OUTER_RADIUS, "Outermost ending radius of chromosome track", DEFAULT_OUTER_RADIUS);
        configBuilder.addDecimal(INNER_RADIUS, "Innermost starting radius of minor-allele ploidy track", DEFAULT_INNER_RADIUS);
        configBuilder.addDecimal(GAP_RADIUS, "Radial gap between tracks", DEFAULT_GAP_RADIUS);
        configBuilder.addDecimal(EXON_RANK_RADIUS, "Radial gap left for exon rank labels", DEFAULT_EXON_RANK_RADIUS);

        configBuilder.addInteger(MAX_PLOT_SVS, "Max SVs on a single plot (set to zero to show unlimited)", DEFAULT_MAX_PLOT_SV_COUNT);

        configBuilder.addInteger(FUSION_HEIGHT, "Height of each fusion in pixels", DEFAULT_FUSION_HEIGHT);
        configBuilder.addInteger(FUSION_LEGEND_ROWS, "Number of rows in protein domain legend ", DEFAULT_FUSION_LEGEND_ROWS);

        configBuilder.addInteger(
                FUSION_LEGEND_HEIGHT_PER_ROW,
                "Height of each row in protein domain legend", DEFAULT_FUSION_LEGEND_HEIGHT_PER_ROW);

        configBuilder.addInteger(CHR_RANGE_HEIGHT, "Chromosome range row height in pixels", DEFAULT_CHR_RANGE_HEIGHT);
        configBuilder.addInteger(CHR_RANGE_COLUMNS, "Chromosome range row columns", DEFAULT_CHR_RANGE_COLUMNS);

        configBuilder.addDecimal(
                GENE_RELATIVE_SIZE,
                                "Size of gene track relative to segments and copy number alterations", DEFAULT_GENE_RELATIVE_SIZE);
        configBuilder.addDecimal(
                SEGMENT_RELATIVE_SIZE,"Size of segment track relative to copy number alterations and genes",
                DEFAULT_SEGMENT_RELATIVE_SIZE);

        configBuilder.addDecimal(
                CNA_RELATIVE_SIZE, "Size of gene copy number alteration relative to genes and segments", DEFAULT_CNA_RELATIVE_SIZE);

        configBuilder.addInteger(MIN_LABEL_SIZE, "Minimum label size in pixels", DEFAULT_MIN_LABEL_SIZE);
        configBuilder.addInteger(MAX_LABEL_SIZE, "Maximum label size in pixels", DEFAULT_MAX_LABEL_SIZE);
        configBuilder.addInteger(MAX_DISTANCE_LABELS, "Maximum number of distance labels before removing them",
                DEFAULT_MAX_DISTANCE_LABELS);
        configBuilder.addInteger(MAX_POSITION_LABELS,
                "Maximum number of position labels before increasing distance between labels", DEFAULT_MAX_POSITION_LABELS);

        configBuilder.addInteger(
                MAX_GENE_CHARACTERS, "Maximum number of character in gene allowed before scaling", DEFAULT_MAX_GENE_CHARACTERS);

        configBuilder.addInteger(MIN_LINE_SIZE, "Minimum line size in pixels", DEFAULT_MIN_LINE_SIZE);
        configBuilder.addInteger(MAX_LINE_SIZE, "Maximum line size in pixels", DEFAULT_MAX_LINE_SIZE);
        configBuilder.addInteger(GLYPH_SIZE, "Size of glyphs in pixels", DEFAULT_GLYPH_SIZE);

        configBuilder.addFlag(INTERPOLATE_CNA_POSITIONS, "Interpolate copy number positions rather than adjust scale");
        configBuilder.addFlag(INTERPOLATE_EXON_POSITIONS, "Interpolate exon positions rather than adjust scale");

        configBuilder.addFlag(EXACT_POSITION, "Display exact positions at all break ends");
        configBuilder.addFlag(SHOW_SV_ID, "Display SV Id next to position");
        configBuilder.addFlag(STEP, "Create image at each step of derivative chromosome");
    }

    public long labelSize(long count)
    {
        if(count > MaxNumberOfDistanceLabels)
        {
            return MinLabelSize;
        }

        return Math.round(MaxLabelSize - 1d * count * (MaxLabelSize - MinLabelSize) / MaxNumberOfDistanceLabels);
    }

    public CircosConfig(final ConfigBuilder configBuilder)
    {
        MinLabelSize = configBuilder.getInteger(MIN_LABEL_SIZE);
        MaxLabelSize = configBuilder.getInteger(MAX_LABEL_SIZE);

        OuterRadius = configBuilder.getDecimal(OUTER_RADIUS);
        InnerRadius = configBuilder.getDecimal(INNER_RADIUS);
        GapRadius = configBuilder.getDecimal(GAP_RADIUS);
        ExonRankRadius = configBuilder.getDecimal(EXON_RANK_RADIUS);

        GeneRelativeSize = configBuilder.getDecimal(GENE_RELATIVE_SIZE);
        SegmentRelativeSize = configBuilder.getDecimal(SEGMENT_RELATIVE_SIZE);
        CopyNumberRelativeSize = configBuilder.getDecimal(CNA_RELATIVE_SIZE);

        MaxPlotSvCount = configBuilder.getInteger(MAX_PLOT_SVS);

        MaxGeneCharacters = configBuilder.getInteger(MAX_GENE_CHARACTERS);
        MaxNumberOfDistanceLabels = configBuilder.getInteger(MAX_DISTANCE_LABELS);
        MaxNumberOfPositionLabels = configBuilder.getInteger(MAX_POSITION_LABELS);

        FusionLegendRows = configBuilder.getInteger(FUSION_LEGEND_ROWS);
        FusionLegendHeightPerRow = configBuilder.getInteger(FUSION_LEGEND_HEIGHT_PER_ROW);
        FusionHeight = configBuilder.getInteger(FUSION_HEIGHT);

        ChromosomeRangeColumns = configBuilder.getInteger(CHR_RANGE_COLUMNS);
        ChromosomeRangeHeight = configBuilder.getInteger(CHR_RANGE_HEIGHT);

        MinLineSize = configBuilder.getInteger(MIN_LINE_SIZE);
        MaxLineSize = configBuilder.getInteger(MAX_LINE_SIZE);
        GlyphSize = configBuilder.getInteger(GLYPH_SIZE);

        InterpolateCopyNumberPositions = configBuilder.hasFlag(INTERPOLATE_CNA_POSITIONS);
        InterpolateExonPositions = configBuilder.hasFlag(INTERPOLATE_EXON_POSITIONS);

        ExactPosition = configBuilder.hasFlag(EXACT_POSITION);
        ShowSvId = configBuilder.hasFlag(SHOW_SV_ID);
        Step = configBuilder.hasFlag(STEP);
    }

    public boolean isValid()
    {
        if(MinLabelSize <= 0)
        {
            VIS_LOGGER.error("Parameter " + MIN_LABEL_SIZE + " should be > 0");
            return false;
        }

        if(MaxLabelSize < MinLabelSize)
        {
            VIS_LOGGER.error("Parameter " + MAX_LABEL_SIZE + " should be > " + MIN_LABEL_SIZE);
            return false;
        }

        if(MinLineSize <= 0)
        {
            VIS_LOGGER.error("Parameter " + MIN_LINE_SIZE + " should be > 0");
            return false;
        }

        if(MaxLineSize < MinLineSize)
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

    public boolean exceedsMaxPlotSvCount(int svCount) { return MaxPlotSvCount > 0 && svCount > MaxPlotSvCount; }

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
