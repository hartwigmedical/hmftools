package com.hartwig.hmftools.svanalysis.visualisation.data;

import static com.hartwig.hmftools.svanalysis.visualisation.circos.Span.maxPositionPerChromosome;
import static com.hartwig.hmftools.svanalysis.visualisation.circos.Span.minPositionPerChromosome;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.refgenome.RefGenome;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.svanalysis.visualisation.circos.SegmentTerminal;
import com.hartwig.hmftools.svanalysis.visualisation.circos.Span;

import org.jetbrains.annotations.NotNull;

public class Segments {

    private static final String HEADER = "SampleId";
    private static final String COMMENT = "#";
    private static final String DELIMITER = ",";
    private static final RefGenome REF_GENOME = RefGenome.HG19;

    @NotNull
    public static Segment entireChromosome(@NotNull final String sampleId, @NotNull final String chromosome) {
        return ImmutableSegment.builder()
                .sampleId(sampleId)
                .clusterId(-1)
                .chainId(-1)
                .chromosome(chromosome)
                .start(1)
                .end(REF_GENOME.lengths().get(HumanChromosome.fromString(chromosome)))
                .track(0)
                .startTerminal(SegmentTerminal.TELOMERE)
                .endTerminal(SegmentTerminal.TELOMERE)
                .traverseCount(0)
                .build();
    }

    @NotNull
    private static Segment centromere(@NotNull final String sampleId, @NotNull final String chromosome) {
        long position = REF_GENOME.centromeres().get(HumanChromosome.fromString(chromosome));
        return ImmutableSegment.builder()
                .sampleId(sampleId)
                .clusterId(-1)
                .chainId(-1)
                .chromosome(chromosome)
                .start(position)
                .end(position)
                .track(0)
                .startTerminal(SegmentTerminal.CENTROMERE)
                .endTerminal(SegmentTerminal.CENTROMERE)
                .traverseCount(0)
                .build();
    }

    @NotNull
    public static List<Segment> readTracks(@NotNull final String fileName) throws IOException {
        return fromString(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    public static List<Segment> extendTerminals(long terminalDistance, @NotNull final List<Segment> segments,
            @NotNull final List<Link> links, @NotNull final List<GenomePosition> allPositions) {
        final Map<Chromosome, Long> centromeres = REF_GENOME.centromeres();
        final Map<Chromosome, Long> lengths = REF_GENOME.lengths();

        final Map<String, Long> minPositionPerChromosome = minPositionPerChromosome(allPositions);
        final Map<String, Long> maxPositionPerChromosome = maxPositionPerChromosome(allPositions);

        final List<Segment> result = Lists.newArrayList();
        for (Segment segment : segments) {
            final long centromere = centromeres.get(HumanChromosome.fromString(segment.chromosome()));
            final long length = lengths.get(HumanChromosome.fromString(segment.chromosome()));

            if (segment.startTerminal() != SegmentTerminal.NONE) {
                final long minPositionOnChromosome = minPositionPerChromosome.get(segment.chromosome());
                final long startPosition = segment.startTerminal() == SegmentTerminal.CENTROMERE && minPositionOnChromosome < centromere
                        ? centromere
                        : Math.max(1, minPositionOnChromosome - terminalDistance);

                segment = ImmutableSegment.builder().from(segment).start(startPosition).build();
            }

            if (segment.endTerminal() != SegmentTerminal.NONE) {
                final long maxPositionOnChromosome = maxPositionPerChromosome.get(segment.chromosome());
                final long endPosition = segment.endTerminal() == SegmentTerminal.CENTROMERE && maxPositionOnChromosome > centromere
                        ? centromere
                        : Math.min(length, maxPositionOnChromosome + terminalDistance);

                segment = ImmutableSegment.builder().from(segment).end(endPosition).build();
            }

            result.add(segment);
        }

        return incrementOnChromosome(addCentromeres(result), links);
    }

    @NotNull
    public static List<Segment> addCentromeres(@NotNull final List<Segment> segments) {
        if (segments.isEmpty()) {
            return segments;
        }
        final List<Segment> result = Lists.newArrayList(segments);
        final Set<String> existingCentromeres = segments.stream()
                .filter(x -> x.startTerminal() == SegmentTerminal.CENTROMERE || x.endTerminal() == SegmentTerminal.CENTROMERE)
                .map(GenomeRegion::chromosome)
                .collect(Collectors.toSet());

        final Set<String> requiredCentomeres = Sets.newHashSet();
        final List<GenomeRegion> segmentSpan = Span.spanRegions(segments);
        for (final GenomeRegion genomeRegion : segmentSpan) {
            long centromere = REF_GENOME.centromeres().get(HumanChromosome.fromString(genomeRegion.chromosome()));
            if (genomeRegion.start() < centromere && genomeRegion.end() > centromere) {
                requiredCentomeres.add(genomeRegion.chromosome());
            }
        }

        requiredCentomeres.removeAll(existingCentromeres);
        if (!requiredCentomeres.isEmpty()) {
            final String sampleId = segments.get(0).sampleId();
            for (final String requiredCentromere : requiredCentomeres) {
                result.add(centromere(sampleId, requiredCentromere));
            }
        }

        return result;
    }

    @VisibleForTesting
    @NotNull
    static List<Segment> fromString(@NotNull final List<String> lines) {
        final List<Segment> result = Lists.newArrayList();

        for (final String line : lines) {

            if (!line.startsWith(COMMENT) && !line.startsWith(HEADER)) {
                String[] values = line.split(DELIMITER);
                final String start = values[4];
                final String end = values[5];

                Segment newSegment = ImmutableSegment.builder()
                        .sampleId(values[0])
                        .clusterId(Integer.valueOf(values[1]))
                        .chainId(Integer.valueOf(values[2]))
                        .chromosome(values[3])
                        .start(SegmentTerminal.fromString(start) == SegmentTerminal.NONE ? Long.valueOf(start) : Long.valueOf(end))
                        .end(SegmentTerminal.fromString(end) == SegmentTerminal.NONE ? Long.valueOf(end) : Long.valueOf(start))
                        .track(0)
                        .startTerminal(SegmentTerminal.fromString(start))
                        .endTerminal(SegmentTerminal.fromString(end))
                        .traverseCount(Integer.valueOf(values[6]))
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

    @NotNull
    private static List<Segment> incrementOnChromosome3(@NotNull final List<Segment> segments, @NotNull final List<Link> links) {

        final Set<Integer> clustersWithoutSegments =
                links.stream().filter(Link::connectorsOnly).map(Link::clusterId).collect(Collectors.toSet());

        final Map<String, Integer> trackMap = Maps.newHashMap();
        final List<Segment> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final Segment segment : segments) {
            if (segment.clusterId() == -1) {
                result.add(ImmutableSegment.builder().from(segment).track(0).build());
            } else if (!clustersWithoutSegments.contains(segment.clusterId())) {
                final String chromosome = segment.chromosome();
                if (!trackMap.containsKey(chromosome)) {
                    trackMap.put(chromosome, currentTrack);
                } else {
                    currentTrack = Math.max(currentTrack, trackMap.get(chromosome) + 1);
                    trackMap.put(chromosome, currentTrack);
                }

                result.add(ImmutableSegment.builder().from(segment).track(currentTrack).build());
            }
        }

        return result;

    }

    @NotNull
    private static List<Segment> incrementOnChromosome(@NotNull final List<Segment> segments, @NotNull final List<Link> links) {

        final Set<Integer> clustersWithoutSegments =
                links.stream().filter(Link::connectorsOnly).map(Link::clusterId).collect(Collectors.toSet());

        final Map<String, Integer> trackMap = Maps.newHashMap();
        final List<Segment> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final Segment segment : segments) {
            if (segment.clusterId() == -1) {
                result.add(ImmutableSegment.builder().from(segment).track(0).build());
            } else if (!clustersWithoutSegments.contains(segment.clusterId())) {
                final String chromosome = segment.chromosome();
                if (!trackMap.containsKey(chromosome)) {
                    currentTrack = 1;
                } else {
                    currentTrack = trackMap.get(chromosome) + 1;
                }
                trackMap.put(chromosome, currentTrack);

                result.add(ImmutableSegment.builder().from(segment).track(currentTrack).build());
            }
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

}
