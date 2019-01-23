package com.hartwig.hmftools.svanalysis.visualisation;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ScalePosition {

    private static final Logger LOGGER = LogManager.getLogger(ScalePosition.class);

    private final Map<String, Map<Long, Integer>> chromosomePositionMap = Maps.newHashMap();

    public ScalePosition(@NotNull final List<? extends GenomeRegion> regions) {
        this(1, regions);
    }

    public ScalePosition(final int start, @NotNull final List<? extends GenomeRegion> regions) {
        final Set<String> contigs = regions.stream().map(GenomeRegion::chromosome).collect(Collectors.toSet());
        for (final String contig : contigs) {
            final List<Long> contigPositions = Lists.newArrayList();
            regions.stream().filter(x -> x.chromosome().equals(contig)).forEach(x -> {
                contigPositions.add(x.start());
                contigPositions.add(x.end());
            });

            chromosomePositionMap.put(contig, positionMap(start, contigPositions));
        }
    }

    @NotNull
    public List<GenomePosition> original() {
        final List<GenomePosition> result = Lists.newArrayList();

        for (String contig : chromosomePositionMap.keySet()) {
            for (Long position : chromosomePositionMap.get(contig).keySet()) {
                result.add(GenomePositions.create(contig, position));
            }
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    public List<GenomePosition> scaled() {
        final List<GenomePosition> result = Lists.newArrayList();

        for (String contig : chromosomePositionMap.keySet()) {
            for (Integer position : chromosomePositionMap.get(contig).values()) {
                result.add(GenomePositions.create(contig, position));
            }
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    public List<Track> scaleTracks(@NotNull final List<Track> tracks) {
        return tracks.stream().map(x -> scale(x, chromosomePositionMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    public List<Link> scaleLinks(@NotNull final List<Link> links) {
        final List<Link> results = Lists.newArrayList();

        for (final Link link : links) {

            try {
                final ImmutableLink.Builder builder = ImmutableLink.builder().from(link);
                if (HumanChromosome.contains(link.startChromosome())) {
                    builder.startPosition(chromosomePositionMap.get(link.startChromosome()).get(link.startPosition()));
                }

                if (HumanChromosome.contains(link.endChromosome())) {
                    builder.endPosition(chromosomePositionMap.get(link.endChromosome()).get(link.endPosition()));
                }

                results.add(builder.build());
            } catch (Exception e) {
                LOGGER.error("Unable to scale link {}", link);
                throw e;
            }
        }

        return results;
    }

    @NotNull
    public List<CopyNumberAlteration> scaleAlterations(@NotNull final List<CopyNumberAlteration> links) {
        return links.stream().map(x -> scale(x, chromosomePositionMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    private static CopyNumberAlteration scale(@NotNull final CopyNumberAlteration victim, @NotNull final Map<Long, Integer> positionMap) {
        return ImmutableCopyNumberAlteration.builder()
                .from(victim)
                .start(positionMap.get(victim.start()))
                .end(positionMap.get(victim.end()))
                .build();
    }

    @NotNull
    private static Track scale(@NotNull final Track victim, @NotNull final Map<Long, Integer> positionMap) {
        return ImmutableTrack.builder().from(victim).start(positionMap.get(victim.start())).end(positionMap.get(victim.end())).build();
    }

    @VisibleForTesting
    static Map<Long, Integer> positionMap(int start, @NotNull final Long... positionArray) {
        return positionMap(start, Lists.newArrayList(positionArray));
    }

    @NotNull
    private static Map<Long, Integer> positionMap(int start, @NotNull final List<Long> positions) {
        final Map<Long, Integer> results = Maps.newHashMap();
        final List<Long> sortedDistinctPositions = positions.stream().sorted().distinct().collect(Collectors.toList());

        if (!sortedDistinctPositions.isEmpty()) {
            int logPosition = start;
            long lastPosition = sortedDistinctPositions.get(0);
            results.put(lastPosition, logPosition);

            for (int i = 1; i < sortedDistinctPositions.size(); i++) {
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

    static int logDistance(long distance) {
        return (int) Math.floor(Math.pow(Math.log10(distance), 3)) + 1;
    }

}
