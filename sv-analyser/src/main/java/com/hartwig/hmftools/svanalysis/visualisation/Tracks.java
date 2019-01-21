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

    public static List<Track> addLinkTerminals(long distance, @NotNull final List<Track> tracks, @NotNull final List<Link> links) {
        final List<Track> result = Lists.newArrayList();

        final List<Integer> chainIds = links.stream().map(Link::chainId).distinct().collect(Collectors.toList());
        for (Integer chainId : chainIds) {
            final List<Link> chainLinks = links.stream().filter(x -> x.chainId() == chainId).collect(Collectors.toList());
            final List<Track> chainTracks = tracks.stream().filter(x -> x.chainId() == chainId).collect(Collectors.toList());
            result.addAll(addLinkTerminals(chainId, distance, chainTracks, chainLinks));
        }

        return TRACK_INCREMENTER.apply(result);
    }

    @NotNull
    private static List<Track> addLinkTerminals(int chainId, long distance, @NotNull final List<Track> tracks,
            @NotNull final List<Link> links) {
        final List<Track> result = Lists.newArrayList(tracks);

        if (!links.isEmpty()) {
            final Link firstLink = links.get(0);
            if (HumanChromosome.contains(firstLink.startChromosome()) && firstLink.startPosition() != -1) {
                final String contig = firstLink.startChromosome();
                final Optional<Track> matchingTrack = findTrack(GenomePositions.create(contig, firstLink.startPosition()), tracks);
                if (!matchingTrack.isPresent()) {
                    long start = firstLink.startPosition() - firstLink.startOrientation() * distance;
                    result.add(0, create(chainId, contig, start, firstLink.startPosition()));

                }
            }

            final Link finalLink = links.get(links.size() - 1);
            if (HumanChromosome.contains(finalLink.endChromosome()) && finalLink.endPosition() != -1) {
                final String contig = finalLink.endChromosome();
                final Optional<Track> matchingTrack = findTrack(GenomePositions.create(contig, finalLink.endPosition()), tracks);
                if (!matchingTrack.isPresent()) {
                    long end = finalLink.endPosition() - finalLink.endOrientation() * distance;
                    result.add(create(chainId, contig, finalLink.endPosition(), end));

                }
            }

        }

        return result;
    }

    @NotNull
    private static Track create(int chainId, String contig, long start, long end) {
        return ImmutableTrack.builder()
                .chainId(chainId)
                .chromosome(contig)
                .start(Math.min(start, end))
                .end(Math.max(start, end))
                .track(0)
                .build();
    }

    @VisibleForTesting
    @NotNull
    static List<Track> fromString(@NotNull final List<String> lines) {
        final List<Track> result = Lists.newArrayList();

        for (final String line : lines) {

            if (!line.startsWith(COMMENT)) {
                String[] values = line.split(DELIMITER);
                final int chainId = Integer.valueOf(values[0]);
                final String chromosome = values[1];
                final long start = Long.valueOf(values[2]);
                final long end = Long.valueOf(values[3]);

                result.add(create(chainId, chromosome, start, end));

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

}
