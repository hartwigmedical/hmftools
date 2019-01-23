package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class Tracks {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";
    private static final Function<List<Track>, List<Track>> TRACK_INCREMENTER = Tracks::incrementOnChromosome;

    @NotNull
    public static List<Track> readTracksFromFile(@NotNull final String fileName) throws IOException {
        return incrementOnChromosome(fromString(Files.readAllLines(new File(fileName).toPath())));
    }

    @NotNull
    public static Optional<Track> findTrackStart(@NotNull final GenomePosition position, @NotNull final List<Track> tracks) {
        return tracks.stream().filter(x -> x.chromosome().equals(position.chromosome()) && x.start() == position.position()).findFirst();
    }

    @NotNull
    public static Optional<Track> findTrackEnd(@NotNull final GenomePosition position, @NotNull final List<Track> tracks) {
        return tracks.stream().filter(x -> x.chromosome().equals(position.chromosome()) && x.end() == position.position()).findFirst();
    }

    @NotNull
    public static Optional<Track> findTrack(@NotNull final GenomePosition position, @NotNull final List<Track> tracks) {
        final Optional<Track> result = findTrackStart(position, tracks);
        return result.isPresent() ? result : findTrackEnd(position, tracks);
    }

    public static long tracksConnectedTo(int minTrackValue, @NotNull final GenomePosition position, @NotNull final List<Track> tracks) {
        return tracks.stream()
                .filter(x -> x.chromosome().equals(position.chromosome()) && (x.start() == position.position()
                        || x.end() == position.position()) && x.track() >= minTrackValue)
                .count();
    }

    public static List<Track> addMissingTracks(long distance, @NotNull final List<Track> tracks,  @NotNull final List<Link> links) {
        final List<Track> result = Lists.newArrayList();

        final List<Integer> chainIds = links.stream().map(Link::chainId).distinct().collect(Collectors.toList());
        for (Integer chainId : chainIds) {
            final List<Link> chainLinks = links.stream().filter(x -> x.chainId() == chainId).collect(Collectors.toList());
            final List<Track> chainTracks = tracks.stream().filter(x -> x.chainId() == chainId).collect(Collectors.toList());
            result.addAll(addMissingTracks(chainId, distance, chainTracks, chainLinks));
        }

        return TRACK_INCREMENTER.apply(result);
    }

    @NotNull
    private static List<Track> addMissingTracks(int chainId, long distance, @NotNull final List<Track> tracks,
            @NotNull final List<Link> links) {
        final List<Track> result = Lists.newArrayList(tracks);

        int i = 0;

        for (final Link link : links) {
            final String startContig = link.startChromosome();
            final String endContig = link.endChromosome();
            if (HumanChromosome.contains(startContig) && HumanChromosome.contains(endContig)) {
                final GenomePosition linkStart = GenomePositions.create(startContig, link.startPosition());
                final GenomePosition linkEnd = GenomePositions.create(endContig, link.endPosition());
                long tracksConnectedToStart = tracksConnectedTo(0, linkStart, tracks);
                long tracksConnectedToEnd = tracksConnectedTo(0, linkEnd, tracks);

                if (tracksConnectedToStart == 0 || tracksConnectedToEnd > tracksConnectedToStart) {

                    long minChromosomePosition = minPosition(linkStart, tracks);
                    long maxChromosomePosition = maxPosition(linkStart, tracks);

                    long start = link.startOrientation() > 0 ? minChromosomePosition - distance : maxChromosomePosition + distance;
                    final Track additionalTrack = create(link.clusterId(),
                            chainId,
                            startContig,
                            start,
                            link.startPosition(),
                            link.startOrientation() > 0,
                            link.startOrientation() < 0);
                    result.add(i++, additionalTrack);
                }

                if (tracksConnectedToEnd == 0 || tracksConnectedToStart > tracksConnectedToEnd) {

                    long minChromosomePosition = minPosition(linkEnd, tracks);
                    long maxChromosomePosition = maxPosition(linkEnd, tracks);

                    long end = link.endOrientation() > 0 ? minChromosomePosition - distance : maxChromosomePosition + distance;
                    final Track additionalTrack = create(link.clusterId(),
                            chainId,
                            endContig,
                            link.endPosition(),
                            end,
                            link.endOrientation() > 0,
                            link.endOrientation() < 0);
                    result.add(additionalTrack);
                }

            }

        }

        return result;
    }

    private static long minPosition(@NotNull final GenomePosition position, @NotNull final List<Track> tracks) {
        return tracks.stream()
                .filter(x -> x.chromosome().equals(position.chromosome()))
                .mapToLong(GenomeRegion::start)
                .min()
                .orElse(position.position());
    }

    private static long maxPosition(@NotNull final GenomePosition position, @NotNull final List<Track> tracks) {
        return tracks.stream()
                .filter(x -> x.chromosome().equals(position.chromosome()))
                .mapToLong(GenomeRegion::end)
                .min()
                .orElse(position.position());
    }

    @VisibleForTesting
    @NotNull
    static List<Track> fromString(@NotNull final List<String> lines) {
        final List<Track> result = Lists.newArrayList();

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
    static List<Track> incrementOnChromosome(@NotNull final List<Track> tracks) {
        final Map<String, Integer> trackMap = Maps.newHashMap();
        final List<Track> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final Track track : tracks) {
            final String chromosome = track.chromosome();
            if (!trackMap.containsKey(chromosome)) {
                trackMap.put(chromosome, currentTrack);
            } else {
                currentTrack = Math.max(currentTrack, trackMap.get(chromosome) + 1);
                trackMap.put(chromosome, currentTrack);
            }

            result.add(ImmutableTrack.builder().from(track).track(currentTrack).build());

        }

        return result;

    }

    @VisibleForTesting
    @NotNull
    static List<Track> alwaysIncrement(@NotNull final List<Track> tracks) {
        final List<Track> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final Track track : tracks) {
            result.add(ImmutableTrack.builder().from(track).track(currentTrack++).build());
        }

        return result;
    }

    @NotNull
    private static Track create(int clusterId, int chainId, String contig, long start, long end, boolean openStart, boolean openEnd) {
        return ImmutableTrack.builder()
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

}
