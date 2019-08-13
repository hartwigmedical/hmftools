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
    String GENE_RELATIVE_SIZE = "gene_proportion";

    double DEFAULT_GENE_PROPORTION = 0.1;

    boolean displayPosition();

    double geneRelativeSize();

    double segmentRelativeSize();

    double copyNumberRelativeSize();

    default double gapSize()
    {
        return 0.025;
    }

    default double innerRadius()
    {
        return 0.175;
    }

    default int minDistanceLabelSize()
    {
        return 40;
    }

    default int maxDistanceLabelSize()
    {
        return 60;
    }

    default int maxDistanceLabelCount()
    {
        return 200;
    }

    default long distanceLabelSize(long count)
    {
        if (count > maxDistanceLabelCount())
        {
            return 0;
        }

        return Math.round(maxDistanceLabelSize() - 1d * count * (maxDistanceLabelSize() - minDistanceLabelSize()) / maxDistanceLabelCount());
    }

    static void addOptions(@NotNull Options options)
    {
        options.addOption(GENE_RELATIVE_SIZE, true, "Size of gene track relative to copy number, ploidy and segments");
    }

    @NotNull
    static SvCircosConfig createConfig(@NotNull final CommandLine cmd) throws ParseException
    {
        return ImmutableSvCircosConfig.builder()
                .displayPosition(false)
                .geneRelativeSize(1)
                .segmentRelativeSize(1)
                .copyNumberRelativeSize(2)
                .build();
    }

}
