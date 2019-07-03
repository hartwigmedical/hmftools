package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionBuilder;
import com.hartwig.hmftools.common.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExon;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableCopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableExon;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableFusedExon;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableGene;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableLink;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableSegment;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class ScalePosition
{

    private static final Logger LOGGER = LogManager.getLogger(ScalePosition.class);

    private final Map<String, Map<Long, Integer>> chromosomePositionMap = Maps.newHashMap();

    ScalePosition(@NotNull final List<? extends GenomePosition> regions)
    {
        this(1, regions);
    }

    private ScalePosition(final int start, @NotNull final List<? extends GenomePosition> positions)
    {
        final Set<String> contigs = positions.stream().map(GenomePosition::chromosome).collect(Collectors.toSet());
        for (final String contig : contigs)
        {
            final List<Long> contigPositions = positions.stream()
                    .filter(x -> x.chromosome().equals(contig))
                    .map(GenomePosition::position)
                    .collect(Collectors.toList());
            chromosomePositionMap.put(contig, positionMap(start, contigPositions));
        }
    }

    @NotNull
    public List<GenomePosition> scaled()
    {
        final List<GenomePosition> result = Lists.newArrayList();

        for (String contig : chromosomePositionMap.keySet())
        {
            for (Integer position : chromosomePositionMap.get(contig).values())
            {
                result.add(GenomePositions.create(contig, position));
            }
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    public List<Segment> scaleSegments(@NotNull final List<Segment> segments)
    {
        return scale(segments, x -> ImmutableSegment.builder().from(x));
    }

    @NotNull
    public List<CopyNumberAlteration> scaleAlterations(@NotNull final List<CopyNumberAlteration> segments)
    {
        return scale(segments, x -> ImmutableCopyNumberAlteration.builder().from(x));
    }

    @NotNull
    public List<ProteinDomain> interpolateProteinDomains(@NotNull final List<ProteinDomain> exons)
    {
        return exons.stream().map(x -> interpolate(x, y -> ImmutableProteinDomain.builder().from(x))).collect(Collectors.toList());
    }

    @NotNull
    public List<Exon> interpolateExons(@NotNull final List<Exon> exons)
    {
        return exons.stream().map(x -> interpolate(x, y -> ImmutableExon.builder().from(x))).collect(Collectors.toList());
    }

    @NotNull
    public List<Gene> scaleGene(@NotNull final List<Gene> genes)
    {
        return genes.stream().map(x ->
        {
            Map<Long, Integer> positionMap = chromosomePositionMap.get(x.chromosome());
            return scale(x, y -> ImmutableGene.builder().from(y).namePosition(positionMap.get(y.namePosition())), positionMap);
        }).collect(Collectors.toList());
    }

    public List<FusedExon> scaleFusedExon(@NotNull final List<FusedExon> exons)
    {
        return exons.stream().map(x ->
        {
            Map<Long, Integer> positionMap = chromosomePositionMap.get(x.chromosome());
            return scale(x, y -> ImmutableFusedExon.builder()
                    .from(y)
                    .geneStart(positionMap.get(y.geneStart()))
                    .geneEnd(positionMap.get(y.geneEnd())), positionMap);
        }).collect(Collectors.toList());
    }

    @NotNull
    private <T extends GenomeRegion> T interpolate(@NotNull final T exon, Function<T, GenomeRegionBuilder<T>> builderFunction)
    {
        final Map<Long, Integer> positionMap = chromosomePositionMap.get(exon.chromosome());
        assert (positionMap != null && !positionMap.isEmpty());

        return builderFunction.apply(exon)
                .start(interpolate(exon.start(), positionMap))
                .end(interpolate(exon.end(), positionMap))
                .build();
    }

    static int interpolate(long value, Map<Long, Integer> positionMap)
    {

        final Set<Long> keySet = positionMap.keySet();

        if (positionMap.containsKey(value))
        {
            return positionMap.get(value);
        }

        long minValue = keySet.stream().mapToLong(x -> x).min().orElse(0);
        long maxValue = keySet.stream().mapToLong(x -> x).max().orElse(0);

        long closestToStart = keySet.stream().filter(x -> x < value).mapToLong(x -> x).max().orElse(minValue);
        long closestToEnd = keySet.stream().filter(x -> x > value).mapToLong(x -> x).min().orElse(maxValue);
        if (closestToStart == closestToEnd)
        {
            return positionMap.get(closestToStart);
        }

        double longDistanceProportion = Math.abs(value - closestToStart) / ((double) Math.abs(closestToEnd - closestToStart));

        int clostestIntToStart = positionMap.get(closestToStart);
        int clostestIntToEnd = positionMap.get(closestToEnd);

        return clostestIntToStart + (int) Math.floor(longDistanceProportion * Math.abs(clostestIntToEnd - clostestIntToStart));
    }

    public List<GenomeRegion> scaleRegions(@NotNull final List<GenomeRegion> regions)
    {
        return regions.stream().map(x -> scale(x, chromosomePositionMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    public List<Link> scaleLinks(@NotNull final List<Link> links)
    {
        final List<Link> results = Lists.newArrayList();
        for (final Link link : links)
        {
            try
            {
                final ImmutableLink.Builder builder = ImmutableLink.builder().from(link);
                if (link.isValidStart())
                {
                    builder.startPosition(chromosomePositionMap.get(link.startChromosome()).get(link.startPosition()));
                }

                if (link.isValidEnd())
                {
                    builder.endPosition(chromosomePositionMap.get(link.endChromosome()).get(link.endPosition()));
                }

                results.add(builder.build());
            } catch (Exception e)
            {
                LOGGER.error("Unable to scale link {}", link);
                throw e;
            }
        }

        return results;
    }

    @NotNull
    private static GenomeRegion scale(@NotNull final GenomeRegion region, @NotNull final Map<Long, Integer> positionMap)
    {
        return GenomeRegions.create(region.chromosome(), positionMap.get(region.start()), positionMap.get(region.end()));
    }

    @VisibleForTesting
    static Map<Long, Integer> positionMap(int start, @NotNull final Long... positionArray)
    {
        return positionMap(start, Lists.newArrayList(positionArray));
    }

    @NotNull
    private static Map<Long, Integer> positionMap(int start, @NotNull final List<Long> positions)
    {
        final Map<Long, Integer> results = Maps.newHashMap();
        final List<Long> sortedDistinctPositions = positions.stream().sorted().distinct().collect(Collectors.toList());

        if (!sortedDistinctPositions.isEmpty())
        {
            int logPosition = start;
            long lastPosition = sortedDistinctPositions.get(0);
            results.put(lastPosition, logPosition);

            for (int i = 1; i < sortedDistinctPositions.size(); i++)
            {
                long position = sortedDistinctPositions.get(i);
                long linearDistance = position - lastPosition;
                int logDistance = logDistance(linearDistance);
                logPosition = logPosition + logDistance;
                lastPosition = position;

                results.put(lastPosition, logPosition);
            }
        }

        return results;
    }

    static int logDistance(long distance)
    {
        return (int) Math.floor(Math.pow(Math.log10(distance), 3)) + 10;
    }

    @NotNull
    private <T extends GenomeRegion> List<T> scale(@NotNull final List<T> inputs, Function<T, GenomeRegionBuilder<T>> builderFunction)
    {
        return inputs.stream().map(x -> scale(x, builderFunction, chromosomePositionMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    private static <T extends GenomeRegion> T scale(@NotNull final T victim, Function<T, GenomeRegionBuilder<T>> builderFunction,
            @NotNull final Map<Long, Integer> positionMap)
    {
        return builderFunction.apply(victim)
                .start(positionMap.get(victim.start()))
                .end(positionMap.get(victim.end())).build();
    }
}
