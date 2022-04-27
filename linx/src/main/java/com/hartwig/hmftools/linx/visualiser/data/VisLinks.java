package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

public class VisLinks
{
    public static Optional<VisSvData> findStartLink(final GenomePosition position, List<VisSvData> links)
    {
        return links.stream()
                .filter(x -> x.ChrStart.equals(position.chromosome()) && x.PosStart == position.position())
                .findFirst();
    }

    public static Optional<VisSvData> findEndLink(final GenomePosition position, List<VisSvData> links)
    {
        return links.stream()
                .filter(x -> x.ChrEnd.equals(position.chromosome()) && x.PosEnd == position.position())
                .findFirst();
    }

    public static Optional<VisSvData> findLink(final GenomePosition position, final List<VisSvData> links)
    {
        for(final VisSvData link : links)
        {
            if (link.ChrStart.equals(position.chromosome()) && link.PosStart == position.position())
            {
                return Optional.of(link);
            }

            if (link.ChrEnd.equals(position.chromosome()) && link.PosEnd == position.position())
            {
                return Optional.of(link);
            }
        }

        return Optional.empty();
    }

    public static List<VisSvData> readSvData(
            final String fileName, final String sampleId, final List<ChrBaseRegion> restrictedRegions) throws IOException
    {
        List<VisSvData> fileSVs = VisSvData.read(fileName);

        List<VisSvData> visSvData = Lists.newArrayList();

        for(VisSvData sv : fileSVs)
        {
            if(!sv.SampleId.equals(sampleId))
                continue;

            if(!restrictedRegions.isEmpty())
            {
                if(restrictedRegions.stream().anyMatch(x -> x.containsPosition(sv.ChrStart, sv.PosStart) || x.containsPosition(sv.ChrEnd, sv.PosEnd)))
                {
                    visSvData.add(sv);
                }
            }
            else
            {
                visSvData.add(sv);
            }
        }

        return visSvData;
    }

    public static List<GenomePosition> allPositions(final List<VisSvData> links)
    {
        final List<GenomePosition> results = Lists.newArrayList();

        for (final VisSvData link : links)
        {
            if (link.isValidStart() && link.PosStart != -1)
            {
                results.add(GenomePositions.create(link.ChrStart, link.PosStart));
            }

            if (link.isValidEnd() && link.PosEnd != -1)
            {
                results.add(GenomePositions.create(link.ChrEnd, link.PosEnd));
            }
        }

        Collections.sort(results);

        return results;
    }

    public static List<VisSvData> addFrame(final List<VisSegment> segments, final List<VisSvData> links)
    {
        final List<VisSvData> results = Lists.newArrayList();

        for(VisSvData link : links)
        {
            int minConnectedFrame = segments.stream().filter(x -> VisLinks.connected(link, x)).mapToInt(x -> x.Frame).min().orElse(0);
            VisSvData newLink = VisSvData.from(link);
            newLink.Frame = minConnectedFrame + 1;
            results.add(newLink);
        }

        return results;
    }

    private static boolean connected(final VisSvData link, final VisSegment segment)
    {
        boolean connectedAtStart = segment.chromosome().equals(link.ChrStart) && (segment.start() == link.PosStart
                || segment.end() == link.PosStart);
        boolean connectedAtEnd = segment.chromosome().equals(link.ChrEnd) && (segment.start() == link.PosEnd
                || segment.end() == link.PosEnd);
        return connectedAtStart || connectedAtEnd;
    }

}
