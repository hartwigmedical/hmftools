package com.hartwig.hmftools.linx.visualiser.data;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class Connectors
{

    private final boolean showSimpleSvSegments;

    public Connectors(final boolean showSimpleSvSegments)
    {
        this.showSimpleSvSegments = showSimpleSvSegments;
    }

    @NotNull
    public List<Connector> createConnectors(@NotNull final List<Segment> segments, @NotNull final List<VisSvData> links)
    {
        final List<Connector> result = Lists.newArrayList();

        for (Segment segment : segments)
        {
            final ImmutableConnector.Builder builder = ImmutableConnector.builder()
                    .chromosome(segment.chromosome())
                    .clusterId(segment.clusterId())
                    .chainId(segment.chainId())
                    .track(segment.track());

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
            final Optional<VisSvData> optionalStartPositionLink = VisLinks.findLink(startPosition, links);

            if (optionalStartPositionLink.isPresent())
            {
                double startLinkPloidy = optionalStartPositionLink.get().jcn();
                double startLinkPloidyBeforeSegment = VisSegments.segmentPloidyBefore(segment.track(), startPosition, segments);

                if (startLinkPloidy > 0)
                {
                    result.add(builder.position(segment.start())
                            .ploidy(Math.max(0, startLinkPloidy - startLinkPloidyBeforeSegment))
                            .frame(optionalStartPositionLink.get().frame())
                            .build());
                }
            }

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
            final Optional<VisSvData> optionalEndPositionLink = VisLinks.findLink(endPosition, links);

            if (optionalEndPositionLink.isPresent())
            {
                double endLinkPloidy = optionalEndPositionLink.get().jcn();
                if (endLinkPloidy > 0)
                {
                    double endLinkPloidyBeforeSegment = VisSegments.segmentPloidyBefore(segment.track(), endPosition, segments);
                    result.add(builder.position(segment.end())
                            .ploidy(Math.max(0, endLinkPloidy - endLinkPloidyBeforeSegment))
                            .frame(optionalEndPositionLink.get().frame())
                            .build());
                }
            }
        }

        links.forEach(x -> result.addAll(create(x)));

        return result;
    }

    @NotNull
    private List<Connector> create(@NotNull final VisSvData link)
    {
        @NotNull
        final List<Connector> result = Lists.newArrayList();

        if (link.connectorsOnly(showSimpleSvSegments))
        {
            final ImmutableConnector.Builder builder = ImmutableConnector.builder()
                    .clusterId(link.clusterId())
                    .chainId(link.chainId())
                    .ploidy(link.jcn())
                    .frame(0)
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
