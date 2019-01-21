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

class TrackFile {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";

    @NotNull
    public static List<Track> readLinks(@NotNull final String fileName) throws IOException {
        return maintainTrackBetweenChromosomes(Files.readAllLines(new File(fileName).toPath()));
    }

    @VisibleForTesting
    @NotNull
    static List<Track> maintainTrackBetweenChromosomes(@NotNull final List<String> lines) {
        final Map<String, Integer> trackMap = Maps.newHashMap();
        final List<Track> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final String line : lines) {

            if (!line.startsWith(COMMENT)) {
                String[] values = line.split(DELIMITER);
                final String chromosome = values[0];
                final long start = Long.valueOf(values[1]);
                final long end = Long.valueOf(values[2]);

                if (!trackMap.containsKey(chromosome)) {
                    trackMap.put(chromosome, currentTrack);
                } else {
                    currentTrack = Math.max(currentTrack, trackMap.get(chromosome) + 1);
                    trackMap.put(chromosome, currentTrack);
                }

                result.add(ImmutableTrack.builder().chromosome(chromosome).start(start).end(end).track(currentTrack).build());

            }

        }

        return result;

    }

    @VisibleForTesting
    @NotNull
    static List<Track> alwaysBumpTrack(@NotNull final List<String> lines) {
        final List<Track> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final String line : lines) {

            if (!line.startsWith(COMMENT)) {
                String[] values = line.split(DELIMITER);
                final String chromosome = values[0];
                final long start = Long.valueOf(values[1]);
                final long end = Long.valueOf(values[2]);

                result.add(ImmutableTrack.builder().chromosome(chromosome).start(start).end(end).track(currentTrack++).build());

            }
        }

        return result;
    }

}
