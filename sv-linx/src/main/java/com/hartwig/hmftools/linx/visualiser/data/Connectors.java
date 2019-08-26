package com.hartwig.hmftools.linx.visualiser.data;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class Connectors
{

    @NotNull
    public static List<Connector> createConnectors(@NotNull final List<Segment> segments, @NotNull final List<Link> links)
    {
        final List<Connector> result = Lists.newArrayList();

        for (Segment segment : segments)
        {
            final ImmutableConnector.Builder builder = ImmutableConnector.builder()
                    .chromosome(segment.chromosome())
                    .clusterId(segment.clusterId())
                    .chainId(segment.chainId())
                    .ploidy(segment.ploidy())
                    .track(segment.track());

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
            double startLinkUsage = Links.linkTraverseCount(startPosition, links);

            if (startLinkUsage > 0)
            {
                result.add(builder.position(segment.start()).build());
            }

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
            double endLinkUsage = Links.linkTraverseCount(endPosition, links);

            if (endLinkUsage > 0)
            {
                result.add(builder.position(segment.end()).build());
            }

        }

        links.forEach(x -> result.addAll(create(x)));

        return result;
    }

    @NotNull
    private static List<Connector> create(@NotNull final Link link)
    {
        @NotNull
        final List<Connector> result = Lists.newArrayList();

        if (link.connectorsOnly())
        {
            final ImmutableConnector.Builder builder = ImmutableConnector.builder()
                    .clusterId(link.clusterId())
                    .chainId(link.chainId())
                    .ploidy(link.ploidy())
                    .track(0);

            if (link.isValidStart())
            {
                final Connector start = builder.chromosome(link.startChromosome()).position(link.startPosition()).build();
                result.add(start);
            }

            if (link.isValidEnd())
            {
                final Connector start = builder.chromosome(link.endChromosome()).position(link.endPosition()).build();
                result.add(start);
            }

        }

        return result;
    }

}
