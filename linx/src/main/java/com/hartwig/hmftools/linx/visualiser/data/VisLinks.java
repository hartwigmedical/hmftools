package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;

public class VisLinks
{
    public static Optional<VisSvData> findStartLink(final GenomePosition position, List<VisSvData> links)
    {
        return links.stream()
                .filter(x -> x.startChromosome().equals(position.chromosome()) && x.startPosition() == position.position())
                .findFirst();
    }

    public static Optional<VisSvData> findEndLink(final GenomePosition position, List<VisSvData> links)
    {
        return links.stream()
                .filter(x -> x.endChromosome().equals(position.chromosome()) && x.endPosition() == position.position())
                .findFirst();
    }

    public static Optional<VisSvData> findLink(final GenomePosition position, final List<VisSvData> links)
    {
        for(final VisSvData link : links)
        {
            if (link.startChromosome().equals(position.chromosome()) && link.startPosition() == position.position())
            {
                return Optional.of(link);
            }

            if (link.endChromosome().equals(position.chromosome()) && link.endPosition() == position.position())
            {
                return Optional.of(link);
            }
        }

        return Optional.empty();
    }

    public static List<VisSvData> readSvData(
            final String fileName, final String sampleId, final List<ChrBaseRegion> restrictedRegions) throws IOException
    {
        List<VisSvDataFile> fileSVs = VisSvDataFile.read(fileName);

        List<VisSvData> visSvData = Lists.newArrayList();

        for(VisSvDataFile sv : fileSVs)
        {
            if(!sv.SampleId.equals(sampleId))
                continue;

            if(!restrictedRegions.isEmpty())
            {
                if(restrictedRegions.stream().anyMatch(x -> x.containsPosition(sv.ChrStart, sv.PosStart) || x.containsPosition(sv.ChrEnd, sv.PosEnd)))
                {
                    visSvData.add(VisLinks.fromFile(sv));
                }
            }
            else
            {
                visSvData.add(VisLinks.fromFile(sv));
            }
        }

        return visSvData;
    }

    private static VisSvData fromFile(final VisSvDataFile file)
    {
        return ImmutableVisSvData.builder()
                .sampleId(file.SampleId)
                .clusterId(file.ClusterId)
                .chainId(file.ChainId)
                .svId(file.SvId)
                .type(file.Type)
                .resolvedType(file.ResolvedType)
                .isSynthetic(file.IsSynthetic)
                .startChromosome(file.ChrStart)
                .startPosition(file.PosStart)
                .startOrientation(file.OrientStart)
                .startInfo(file.InfoStart)
                .endChromosome(file.ChrEnd)
                .endPosition(file.PosEnd)
                .endOrientation(file.OrientEnd)
                .endInfo(file.InfoEnd)
                .jcn(file.JCN)
                .inDoubleMinute(file.InDoubleMinute)
                .frame(0)
                .build();
    }

    public static List<GenomePosition> allPositions(final List<VisSvData> links)
    {
        final List<GenomePosition> results = Lists.newArrayList();

        for (final VisSvData link : links)
        {
            if (link.isValidStart() && link.startPosition() != -1)
            {
                results.add(GenomePositions.create(link.startChromosome(), link.startPosition()));
            }

            if (link.isValidEnd() && link.endPosition() != -1)
            {
                results.add(GenomePositions.create(link.endChromosome(), link.endPosition()));
            }
        }

        Collections.sort(results);

        return results;
    }

    public static List<VisSvData> addFrame(final List<Segment> segments, final List<VisSvData> links)
    {
        final List<VisSvData> result = Lists.newArrayList();

        for(VisSvData link : links)
        {
            int minConnectedFrame = segments.stream().filter(x -> VisLinks.connected(link, x)).mapToInt(Segment::frame).min().orElse(0);
            result.add(ImmutableVisSvData.builder().from(link).frame(minConnectedFrame + 1).build());
        }

        return result;
    }

    private static boolean connected(final VisSvData link, final Segment segment)
    {
        boolean connectedAtStart = segment.chromosome().equals(link.startChromosome()) && (segment.start() == link.startPosition()
                || segment.end() == link.startPosition());
        boolean connectedAtEnd = segment.chromosome().equals(link.endChromosome()) && (segment.start() == link.endPosition()
                || segment.end() == link.endPosition());
        return connectedAtStart || connectedAtEnd;
    }

}
