package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.jetbrains.annotations.NotNull;

public class Highlights
{
    private static final String HEADER = "SampleId";
    private static final String COMMENT = "#";
    private static final String DELIMITER = ",";

    @NotNull
    public static List<GenomeRegion> limitHighlightsToRegions(final List<GenomeRegion> highlights, final List<GenomeRegion> segments)
    {
        final List<GenomeRegion> result = Lists.newArrayList();

        for (final GenomeRegion highlight : highlights)
        {
            final String contig = highlight.chromosome();
            final List<GenomeRegion> chromosomeSegments =
                    segments.stream().filter(x -> x.chromosome().equals(contig)).collect(Collectors.toList());
            if (!chromosomeSegments.isEmpty())
            {
                long minTrackPosition = chromosomeSegments.stream().mapToLong(GenomeRegion::start).min().orElse(0);
                long maxTrackPosition = chromosomeSegments.stream().mapToLong(GenomeRegion::end).max().orElse(0);
                if (highlight.end() >= minTrackPosition && highlight.start() <= maxTrackPosition)
                {

                    result.add(GenomeRegions.create(contig,
                            Math.max(minTrackPosition, highlight.start()),
                            Math.min(maxTrackPosition, highlight.end())));
                }
            }
        }
        return result;
    }

    @NotNull
    public static List<GenomeRegion> fragileSites()
    {
        return fromResource("fragile_sites.csv");
    }

    @NotNull
    public static List<GenomeRegion> lineElements()
    {
        return fromResource("line_elements.csv");
    }

    @NotNull
    private static List<GenomeRegion> fromResource(@NotNull final String resource)
    {
        final InputStream inputStream = Highlights.class.getResourceAsStream("/visualisation/" + resource);
        return fromString(new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static List<GenomeRegion> fromString(@NotNull final List<String> lines)
    {
        final List<GenomeRegion> result = Lists.newArrayList();

        for (final String line : lines)
        {

            if (!line.startsWith(COMMENT) && !line.startsWith(HEADER))
            {
                final String[] values = line.split(DELIMITER);
                result.add(GenomeRegions.create(values[0], Long.parseLong(values[1]), Long.parseLong(values[2])));

            }
        }

        return result;
    }

}
