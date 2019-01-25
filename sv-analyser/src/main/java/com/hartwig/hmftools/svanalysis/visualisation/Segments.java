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
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.refgenome.RefGenome;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class Segments {

    private static final String HEADER = "SampleId";
    private static final String COMMENT = "#";
    private static final String DELIMITER = ",";
    private static final Function<List<Segment>, List<Segment>> TRACK_INCREMENTER = Segments::incrementOnChromosome;
    private static final RefGenome REF_GENOME = RefGenome.HG19;

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

    @NotNull
    public static List<Segment> extendTerminals(long additionalDistance, @NotNull final List<Segment> segments, @NotNull final List<Link> links) {
        final Map<Chromosome, Long> centromeres = REF_GENOME.centromeres();

        final List<GenomePosition> allPositions = Lists.newArrayList();
        allPositions.addAll(Segments.allPositions(segments));
        allPositions.addAll(Links.allPositions(links));

        final Map<String, Long> minPositionPerChromosome = minPositionPerChromosome(allPositions);
        final Map<String, Long> maxPositionPerChromosome = maxPositionPerChromosome(allPositions);

        final List<Segment> result = Lists.newArrayList();
        for (Segment segment : segments) {
            final long centromere = centromeres.get(HumanChromosome.fromString(segment.chromosome()));

            if (segment.startTerminal() != SegmentTerminal.NONE) {
                final long minPositionOnChromosome = minPositionPerChromosome.get(segment.chromosome());
                final long startPosition = segment.startTerminal() == SegmentTerminal.CENTROMERE && minPositionOnChromosome < centromere
                        ? centromere
                        : minPositionOnChromosome - additionalDistance;

                segment = ImmutableSegment.builder()
                        .from(segment)
                        .start(startPosition)
                        .build();
            }

            if (segment.endTerminal() != SegmentTerminal.NONE) {
                final long maxPositionOnChromosome = maxPositionPerChromosome.get(segment.chromosome());
                final long endPosition = segment.endTerminal() == SegmentTerminal.CENTROMERE && maxPositionOnChromosome > centromere
                        ? centromere
                        : maxPositionOnChromosome + additionalDistance;

                segment = ImmutableSegment.builder()
                        .from(segment)
                        .end(endPosition)
                        .build();
            }

            result.add(segment);
        }

        return TRACK_INCREMENTER.apply(result);
    }

    @VisibleForTesting
    @NotNull
    static List<Segment> fromString(@NotNull final List<String> lines) {
        final List<Segment> result = Lists.newArrayList();

        for (final String line : lines) {

            if (!line.startsWith(COMMENT) && !line.startsWith(HEADER)) {
                String[] values = line.split(DELIMITER);
                final int clusterId = Integer.valueOf(values[1]);
                final int chainId = Integer.valueOf(values[2]);
                final String chromosome = values[3];
                final String start = values[4];
                final String end = values[5];
                final int traverseCount = Integer.valueOf(values[6]);

                Segment newSegment = ImmutableSegment.builder()
                        .clusterId(clusterId)
                        .chainId(chainId)
                        .chromosome(chromosome)
                        .start(SegmentTerminal.fromString(start) == SegmentTerminal.NONE ? Long.valueOf(start) : Long.valueOf(end))
                        .end(SegmentTerminal.fromString(end) == SegmentTerminal.NONE ? Long.valueOf(end) : Long.valueOf(start))
                        .track(0)
                        .startTerminal(SegmentTerminal.fromString(start))
                        .endTerminal(SegmentTerminal.fromString(end))
                        .traverseCount(traverseCount)
                        .build();

                result.add(newSegment);

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
    private static Map<String, Long> maxPositionPerChromosome(@NotNull final List<GenomePosition> tracks) {
        return tracks.stream().collect(Collectors.toMap(GenomePosition::chromosome, GenomePosition::position, Math::max));
    }

    @NotNull
    private static Map<String, Long> minPositionPerChromosome(@NotNull final List<GenomePosition> tracks) {
        return tracks.stream().collect(Collectors.toMap(GenomePosition::chromosome, GenomePosition::position, Math::min));
    }

    @NotNull
    public static List<GenomePosition> allPositions(@NotNull final List<? extends GenomeRegion> segments) {
        final List<GenomePosition> results = Lists.newArrayList();

        for (final GenomeRegion segment : segments) {
            if (HumanChromosome.contains(segment.chromosome())) {
                results.add(GenomePositions.create(segment.chromosome(), segment.start()));
                results.add(GenomePositions.create(segment.chromosome(), segment.end()));
            }
        }

        Collections.sort(results);
        return results;
    }

}
