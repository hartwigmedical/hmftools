package com.hartwig.hmftools.svvisualise.data;

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
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.refgenome.RefGenome;
import com.hartwig.hmftools.svvisualise.circos.SegmentTerminal;
import com.hartwig.hmftools.svvisualise.circos.Span;

import org.jetbrains.annotations.NotNull;

public class Segments {

    private static final String HEADER = "SampleId";
    private static final String COMMENT = "#";
    private static final String DELIMITER = ",";
    private static final RefGenome REF_GENOME = RefGenome.HG19;

    @NotNull
    public static Segment chromosome(@NotNull final String sampleId, @NotNull final String chromosome) {
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
    public static List<Segment> readTracks(@NotNull final String fileName) throws IOException {
        return fromString(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    public static List<Segment> ensureCoverage(long terminalDistance, @NotNull final List<Segment> segments,
            @NotNull final List<Link> links, @NotNull final List<Exon> exons) {
        final Map<Chromosome, Long> centromeres = REF_GENOME.centromeres();

        final List<GenomePosition> allPositions = Lists.newArrayList();
        allPositions.addAll(Span.allPositions(segments));
        allPositions.addAll(Links.allPositions(links));
        allPositions.addAll(Span.allPositions(exons));

        final Map<String, Long> minPositionPerChromosome = minPositionPerChromosome(allPositions);
        final Map<String, Long> maxPositionPerChromosome = maxPositionPerChromosome(allPositions);

        final List<Segment> result = Lists.newArrayList();
        for (Segment segment : segments) {
            final long centromere = centromeres.get(HumanChromosome.fromString(segment.chromosome()));

            if (segment.startTerminal() != SegmentTerminal.NONE) {
                final long minPositionOnChromosome = minPositionPerChromosome.get(segment.chromosome());
                final long startPosition = segment.startTerminal() == SegmentTerminal.CENTROMERE && minPositionOnChromosome < centromere
                        ? centromere
                        : minPositionOnChromosome - terminalDistance;

                segment = ImmutableSegment.builder().from(segment).start(startPosition).build();
            }

            if (segment.endTerminal() != SegmentTerminal.NONE) {
                final long maxPositionOnChromosome = maxPositionPerChromosome.get(segment.chromosome());
                final long endPosition = segment.endTerminal() == SegmentTerminal.CENTROMERE && maxPositionOnChromosome > centromere
                        ? centromere
                        : maxPositionOnChromosome + terminalDistance;

                segment = ImmutableSegment.builder().from(segment).end(endPosition).build();
            }

            result.add(segment);
        }

        return incrementOnChromosome(result, links);
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

    @VisibleForTesting
    @NotNull
    static List<Segment> incrementOnChromosome(@NotNull final List<Segment> segments, @NotNull final List<Link> links) {

        final Set<Integer> simpleClusters = links.stream().filter(Link::isSimpleSV).map(Link::clusterId).collect(Collectors.toSet());

        final Map<String, Integer> trackMap = Maps.newHashMap();
        final List<Segment> result = Lists.newArrayList();

        int currentTrack = 1;
        for (final Segment segment : segments) {
            if (simpleClusters.contains(segment.clusterId()) || segment.clusterId() == -1) {
                result.add(ImmutableSegment.builder().from(segment).track(0).build());
            } else {

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
}
