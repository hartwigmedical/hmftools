package com.hartwig.hmftools.linx.visualisation.data;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.linx.visualiser.VisSvDataFile;

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
        final Optional<Link> result = findStartLink(position, links);
        return result.isPresent() ? result : findEndLink(position, links);
    }

    public static int linkTraverseCount(@NotNull final GenomePosition position, @NotNull final List<Link> links)
    {
        return Links.findStartLink(position, links).map(Link::traverseCount).orElse(0) + Links.findEndLink(position, links)
                .map(Link::traverseCount)
                .orElse(0);
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
                .startChromosome(file.ChrStart)
                .startPosition(file.PosStart)
                .startOrientation(file.OrientStart)
                .startInfo(file.InfoStart)
                .endChromosome(file.ChrEnd)
                .endPosition(file.PosEnd)
                .endOrientation(file.OrientEnd)
                .endInfo(file.InfoEnd)
                .traverseCount(file.TraverseCount)
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

}
