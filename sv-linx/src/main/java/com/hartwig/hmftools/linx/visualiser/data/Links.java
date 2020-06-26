package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;

import org.jetbrains.annotations.NotNull;

public class Links
{

    @NotNull
    public static Optional<Link> findStartLink(@NotNull final GenomePosition position, @NotNull List<Link> links)
    {
        return links.stream()
                .filter(x -> x.startChromosome().equals(position.chromosome()) && x.startPosition() == position.position())
                .findFirst();
    }

    @NotNull
    public static Optional<Link> findEndLink(@NotNull final GenomePosition position, @NotNull List<Link> links)
    {
        return links.stream()
                .filter(x -> x.endChromosome().equals(position.chromosome()) && x.endPosition() == position.position())
                .findFirst();
    }

    @NotNull
    public static Optional<Link> findLink(@NotNull final GenomePosition position, @NotNull final List<Link> links)
    {
        for (final Link link : links)
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

    @NotNull
    public static List<Link> clean(@NotNull final List<Link> links)
    {
        return links.stream()
                .filter(x -> HumanChromosome.contains(x.startChromosome()) && HumanChromosome.contains(x.endChromosome()))
                .collect(Collectors.toList());
    }

    @NotNull
    public static List<Link> readLinks(@NotNull final String fileName) throws IOException
    {
        return VisSvDataFile.read(fileName).stream().map(Links::fromFile).collect(Collectors.toList());
    }

    @NotNull
    private static Link fromFile(@NotNull final VisSvDataFile file)
    {
        return ImmutableLink.builder()
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
                .ploidy(file.JCN)
                .inDoubleMinute(file.InDoubleMinute)
                .frame(0)
                .build();
    }

    @NotNull
    public static List<GenomePosition> allPositions(@NotNull final List<Link> links)
    {
        final List<GenomePosition> results = Lists.newArrayList();

        for (final Link link : links)
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

    public static List<Link> addFrame(@NotNull final List<Segment> segments, @NotNull final List<Link> links)
    {
        final List<Link> result = Lists.newArrayList();

        for (Link link : links)
        {
            int minConnectedFrame = segments.stream().filter(x -> Links.connnected(link, x)).mapToInt(Segment::frame).min().orElse(0);
            result.add(ImmutableLink.builder().from(link).frame(minConnectedFrame + 1).build());
        }

        return result;
    }

    private static boolean connnected(@NotNull final Link link, @NotNull final Segment segment)
    {
        boolean connectedAtStart = segment.chromosome().equals(link.startChromosome()) && (segment.start() == link.startPosition()
                || segment.end() == link.startPosition());
        boolean connectedAtEnd = segment.chromosome().equals(link.endChromosome()) && (segment.start() == link.endPosition()
                || segment.end() == link.endPosition());
        return connectedAtStart || connectedAtEnd;
    }

}
