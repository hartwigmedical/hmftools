package com.hartwig.hmftools.linx.visualiser;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SvCircosConfig
{
    String DISPLAY_POSITION = "displayPosition";
    String GENE_RELATIVE_SIZE = "gene_relative_size";
    String CNA_RELATIVE_SIZE = "cna_relative_size";
    String SEGMENT_RELATIVE_SIZE = "segment_relative_size";

    double DEFAULT_GENE_RELATIVE_SIZE = 0.5;
    double DEFAULT_SEGMENT_RELATIVE_SIZE = 1;
    double DEFAULT_CNA_RELATIVE_SIZE = 2;
    double DEFAULT_IDEOGRAM_RADIUS = 0.88;
    double DEFAULT_INNER_RADIUS = 0.2;
    double DEFAULT_GAP_RADIUS = 0.025;

    int MIN_TEXT_SIZE = 35;
    int MAX_TEXT_SIZE = 40;

    boolean displayPosition();

    double geneRelativeSize();

    double segmentRelativeSize();

    double copyNumberRelativeSize();

    default double ideogramRadius()
    {
        return DEFAULT_IDEOGRAM_RADIUS;
    }

    default double gapRadius()
    {
        return DEFAULT_GAP_RADIUS;
    }

    default double innerRadius()
    {
        return DEFAULT_INNER_RADIUS;
    }

    default int minLabelSize()
    {
        return MIN_TEXT_SIZE;
    }

    default int maxLabelSize()
    {
        return MAX_TEXT_SIZE;
    }

    default int maxDistanceLabels()
    {
        return 100;
    }

    default long labelSize(long count)
    {
        if (count > maxDistanceLabels())
        {
            return 0;
        }

        return Math.round(maxLabelSize() - 1d * count * (maxLabelSize() - minLabelSize()) / maxDistanceLabels());
    }

    static void addOptions(@NotNull Options options)
    {
        options.addOption(GENE_RELATIVE_SIZE, true, "Size of gene track relative to copy number, ploidy and segments");
    }

    @NotNull
    static SvCircosConfig createConfig(@NotNull final CommandLine cmd) throws ParseException
    {
        return ImmutableSvCircosConfig.builder()
                .displayPosition(true)
                .geneRelativeSize(DEFAULT_GENE_RELATIVE_SIZE)
                .segmentRelativeSize(DEFAULT_SEGMENT_RELATIVE_SIZE)
                .copyNumberRelativeSize(DEFAULT_CNA_RELATIVE_SIZE)
                .build();
    }

}
