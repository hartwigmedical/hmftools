package com.hartwig.hmftools.linx.visualiser.data;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

import org.jetbrains.annotations.NotNull;

public class Connectors
{
    private final boolean showSimpleSvSegments;

    public Connectors(final boolean showSimpleSvSegments)
    {
        this.showSimpleSvSegments = showSimpleSvSegments;
    }

    @NotNull
    public List<Connector> createConnectors(@NotNull final List<VisSegment> segments, @NotNull final List<VisSvData> links)
    {
        final List<Connector> result = Lists.newArrayList();

        for (VisSegment segment : segments)
        {
            final ImmutableConnector.Builder builder = ImmutableConnector.builder()
                    .chromosome(segment.chromosome())
                    .clusterId(segment.ClusterId)
                    .chainId(segment.ChainId)
                    .track(segment.Track);

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
            final Optional<VisSvData> optionalStartPositionLink = VisLinks.findLink(startPosition, links);

            if (optionalStartPositionLink.isPresent())
            {
                double startLinkPloidy = optionalStartPositionLink.get().JCN;
                double startLinkPloidyBeforeSegment = VisSegments.segmentPloidyBefore(segment.Track, startPosition, segments);

                if (startLinkPloidy > 0)
                {
                    result.add(builder.position(segment.start())
                            .ploidy(Math.max(0, startLinkPloidy - startLinkPloidyBeforeSegment))
                            .frame(optionalStartPositionLink.get().Frame)
                            .build());
                }
            }

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
            final Optional<VisSvData> optionalEndPositionLink = VisLinks.findLink(endPosition, links);

            if (optionalEndPositionLink.isPresent())
            {
                double endLinkPloidy = optionalEndPositionLink.get().JCN;
                if (endLinkPloidy > 0)
                {
                    double endLinkPloidyBeforeSegment = VisSegments.segmentPloidyBefore(segment.Track, endPosition, segments);
                    result.add(builder.position(segment.end())
                            .ploidy(Math.max(0, endLinkPloidy - endLinkPloidyBeforeSegment))
                            .frame(optionalEndPositionLink.get().Frame)
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
                    .clusterId(link.ClusterId)
                    .chainId(link.ChainId)
                    .ploidy(link.JCN)
                    .frame(0)
                    .track(0);

            if (link.isValidStart())
            {
                final Connector start = builder.chromosome(link.ChrStart).position(link.PosStart).build();
                result.add(start);
            }

            if (link.isValidEnd())
            {
                final Connector start = builder.chromosome(link.ChrEnd).position(link.PosEnd).build();
                result.add(start);
            }

        }

        return result;
    }

}
