package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class TrackFactory {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";

    @NotNull
    public static List<Track> readTracks(@NotNull final String fileName) throws IOException {
        return incrementOnChromosome(fromString(Files.readAllLines(new File(fileName).toPath())));
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

                result.add(ImmutableTrack.builder().chainId(chainId).chromosome(chromosome).start(start).end(end).track(0).build());

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
