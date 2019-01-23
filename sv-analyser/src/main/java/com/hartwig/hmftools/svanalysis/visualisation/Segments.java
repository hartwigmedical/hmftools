package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class Segments {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";
    private static final Function<List<Segment>, List<Segment>> TRACK_INCREMENTER = Segments::incrementOnChromosome;

    @NotNull
    public static List<Segment> readTracksFromFile(@NotNull final String fileName) throws IOException {
        return incrementOnChromosome(fromString(Files.readAllLines(new File(fileName).toPath())));
    }

    @NotNull
    public static Optional<Segment> findTrackStart(@NotNull final GenomePosition position, @NotNull final List<Segment> segments) {
        return segments.stream().filter(x -> x.chromosome().equals(position.chromosome()) && x.start() == position.position()).findFirst();
    }

    @NotNull
    public static Optional<Segment> findTrackEnd(@NotNull final GenomePosition position, @NotNull final List<Segment> segments) {
        return segments.stream().filter(x -> x.chromosome().equals(position.chromosome()) && x.end() == position.position()).findFirst();
    }

    @NotNull
    public static Optional<Segment> findTrack(@NotNull final GenomePosition position, @NotNull final List<Segment> segments) {
        final Optional<Segment> result = findTrackStart(position, segments);
        return result.isPresent() ? result : findTrackEnd(position, segments);
    }

    public static long tracksConnectedTo(int minTrackValue, @NotNull final GenomePosition position, @NotNull final List<Segment> segments) {
        return segments.stream()
                .filter(x -> x.chromosome().equals(position.chromosome()) && (x.start() == position.position()
                        || x.end() == position.position()) && x.track() >= minTrackValue)
                .count();
    }

    public static List<Segment> addMissingTracks(long distance, @NotNull final List<Segment> segments, @NotNull final List<Link> links) {

        final List<GenomePosition> allPositions = Lists.newArrayList();
        allPositions.addAll(Segments.allPositions(segments));
        allPositions.addAll(Links.allPositions(links));

        final Map<String, Long> minPositionPerChromosome = minPositionPerChromosome(allPositions);
        final Map<String, Long> maxPositionPerChromosome = maxPositionPerChromosome(allPositions);

        final List<Segment> result = Lists.newArrayList();

        final List<Integer> chainIds = links.stream().map(Link::chainId).distinct().collect(Collectors.toList());
        for (Integer chainId : chainIds) {
            final List<Link> chainLinks = links.stream().filter(x -> x.chainId() == chainId).collect(Collectors.toList());
            final List<Segment> chainSegments = segments.stream().filter(x -> x.chainId() == chainId).collect(Collectors.toList());
            result.addAll(addMissingTracks(chainId, distance, chainSegments, chainLinks, minPositionPerChromosome, maxPositionPerChromosome));
        }

        return TRACK_INCREMENTER.apply(result);
    }

    @NotNull
    private static List<Segment> addMissingTracks(int chainId, long distance, @NotNull final List<Segment> segments,
            @NotNull final List<Link> links, @NotNull final Map<String, Long> minPositionPerChromosome,
            @NotNull final Map<String, Long> maxPositionPerChromosome) {
        final List<Segment> result = Lists.newArrayList(segments);

        int i = 0;

        for (final Link link : links) {
            final String startContig = link.startChromosome();
            final String endContig = link.endChromosome();

            boolean isStartValid = HumanChromosome.contains(startContig);
            boolean isEndValid = HumanChromosome.contains(endContig);

            final GenomePosition linkStart = GenomePositions.create(startContig, link.startPosition());
            final GenomePosition linkEnd = GenomePositions.create(endContig, link.endPosition());

            long tracksConnectedToStart = tracksConnectedTo(0, linkStart, segments);
            long tracksConnectedToEnd = tracksConnectedTo(0, linkEnd, segments);

            if (isStartValid && (tracksConnectedToStart == 0 || tracksConnectedToEnd > tracksConnectedToStart)) {

                long minChromosomePosition = Optional.ofNullable(minPositionPerChromosome.get(startContig)).orElse(link.startPosition());
                long maxChromosomePosition = Optional.ofNullable(maxPositionPerChromosome.get(startContig)).orElse(link.startPosition());

                long start = link.startOrientation() > 0 ? minChromosomePosition - distance : maxChromosomePosition + distance;
                final Segment additionalSegment = create(link.clusterId(),
                        chainId,
                        startContig,
                        start,
                        link.startPosition(),
                        link.startOrientation() > 0,
                        link.startOrientation() < 0);
                result.add(i++, additionalSegment);
            }

            if (isEndValid && (tracksConnectedToEnd == 0 || tracksConnectedToStart > tracksConnectedToEnd)) {

                long minChromosomePosition = Optional.ofNullable(minPositionPerChromosome.get(endContig)).orElse(link.endPosition());
                long maxChromosomePosition = Optional.ofNullable(maxPositionPerChromosome.get(endContig)).orElse(link.endPosition());

                long end = link.endOrientation() > 0 ? minChromosomePosition - distance : maxChromosomePosition + distance;
                final Segment additionalSegment = create(link.clusterId(),
                        chainId,
                        endContig,
                        link.endPosition(),
                        end,
                        link.endOrientation() > 0,
                        link.endOrientation() < 0);
                result.add(additionalSegment);
            }

        }

        return result;
    }

    @VisibleForTesting
    @NotNull
    static List<Segment> fromString(@NotNull final List<String> lines) {
        final List<Segment> result = Lists.newArrayList();

        for (final String line : lines) {

            if (!line.startsWith(COMMENT)) {
                String[] values = line.split(DELIMITER);
                final int clusterId = Integer.valueOf(values[0]);
                final int chainId = Integer.valueOf(values[1]);
                final String chromosome = values[2];
                final long start = Long.valueOf(values[3]);
                final long end = Long.valueOf(values[4]);

                result.add(create(clusterId, chainId, chromosome, start, end, false, false));

            }
        }

        return result;
    }

    @VisibleForTesting
    @NotNull
    static List<Segment> incrementOnChromosome(@NotNull final List<Segment> segments) {
        final Map<String, Integer> trackMap = Maps.newHashMap();
        final List<Segment> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final Segment segment : segments) {
            final String chromosome = segment.chromosome();
            if (!trackMap.containsKey(chromosome)) {
                trackMap.put(chromosome, currentTrack);
            } else {
                currentTrack = Math.max(currentTrack, trackMap.get(chromosome) + 1);
                trackMap.put(chromosome, currentTrack);
            }

            result.add(ImmutableSegment.builder().from(segment).track(currentTrack).build());

        }

        return result;

    }

    @VisibleForTesting
    @NotNull
    static List<Segment> alwaysIncrement(@NotNull final List<Segment> segments) {
        final List<Segment> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final Segment segment : segments) {
            result.add(ImmutableSegment.builder().from(segment).track(currentTrack++).build());
        }

        return result;
    }

    @NotNull
    private static Segment create(int clusterId, int chainId, String contig, long start, long end, boolean openStart, boolean openEnd) {
        return ImmutableSegment.builder()
                .clusterId(clusterId)
                .chainId(chainId)
                .chromosome(contig)
                .start(Math.min(start, end))
                .end(Math.max(start, end))
                .track(0)
                .openStart(openStart)
                .openEnd(openEnd)
                .build();
    }

    @NotNull
    private static Map<String, Long> maxPositionPerChromosome(@NotNull final List<GenomePosition> tracks) {
        return tracks.stream().collect(Collectors.toMap(GenomePosition::chromosome, GenomePosition::position, Math::max));
    }

    @NotNull
    private static Map<String, Long> minPositionPerChromosome(@NotNull final List<GenomePosition> tracks) {
        return tracks.stream().collect(Collectors.toMap(GenomePosition::chromosome, GenomePosition::position, Math::min));
    }

    @NotNull
    public static List<GenomePosition> allPositions(@NotNull final List<Segment> segments) {
        final List<GenomePosition> results = Lists.newArrayList();

        for (final Segment segment : segments) {
            if (HumanChromosome.contains(segment.chromosome())) {
                results.add(GenomePositions.create(segment.chromosome(), segment.start()));
                results.add(GenomePositions.create(segment.chromosome(), segment.end()));
            }

        }

        Collections.sort(results);

        return results;
    }

}
