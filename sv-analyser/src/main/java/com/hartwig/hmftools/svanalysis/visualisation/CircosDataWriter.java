package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class CircosDataWriter {

    private static final String DELIMITER = "\t";

    private final String filePrefix;
    private final int maxTracks;

    public CircosDataWriter(@NotNull final String sample, @NotNull final String outputDir, final int maxTracks) {
        this.filePrefix = outputDir + File.separator + sample;
        this.maxTracks = maxTracks;
    }

    public void write(@NotNull final List<Segment> unadjustedSegments, @NotNull final List<Link> unadjustedLinks,
            @NotNull final List<CopyNumberAlteration> unadjustedAlterations) throws IOException {

        final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
        unadjustedPositions.addAll(Segments.allPositions(unadjustedSegments));
        unadjustedPositions.addAll(Segments.allPositions(unadjustedAlterations));

        final ScalePosition scalePosition = new ScalePosition(unadjustedPositions);
        final List<GenomePosition> scaledPositions = scalePosition.scaled();
        final Map<String, Integer> contigLengths = contigLengths(scaledPositions);

        final List<Segment> segments = scalePosition.scaleTracks(unadjustedSegments);
        final List<Link> links = scalePosition.scaleLinks(unadjustedLinks);
        final List<CopyNumberAlteration> alterations = scalePosition.scaleAlterations(unadjustedAlterations);

        final String textPath = filePrefix + ".text.circos";
        Files.write(new File(textPath).toPath(), createPositionText(scalePosition.original(), scaledPositions));

        final String histogramPath = filePrefix + ".histogram.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(contigLengths, segments));

        final String karyotypePath = filePrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(contigLengths));

        final String connectorPath = filePrefix + ".connector.circos";
        Files.write(new File(connectorPath).toPath(), createConnectors(maxTracks, segments, links));

        final String linkPath = filePrefix + ".link.circos";
        Files.write(new File(linkPath).toPath(), createLinks(segments, links));

        final String scatterPath = filePrefix + ".scatter.circos";
        Files.write(new File(scatterPath).toPath(), createScatter(segments, links));

        final String cnaPath = filePrefix + ".cna.circos";
        Files.write(new File(cnaPath).toPath(), createCNA(alterations));

        final String mapPath = filePrefix + ".map.circos";
        Files.write(new File(mapPath).toPath(), createMinorAllelePloidy(alterations));
    }

    @NotNull
    private List<String> createCNA(@NotNull final List<CopyNumberAlteration> alterations) {
        final List<String> result = Lists.newArrayList();
        for (CopyNumberAlteration alteration : alterations) {
            final String cna = new StringJoiner(DELIMITER).add(circosContig(alteration.chromosome()))
                    .add(String.valueOf(alteration.start()))
                    .add(String.valueOf(alteration.end()))
                    .add(String.valueOf(alteration.copyNumber() - 2))
                    .toString();
            result.add(cna);
        }

        return result;
    }

    @NotNull
    private List<String> createMinorAllelePloidy(@NotNull final List<CopyNumberAlteration> alterations) {
        final List<String> result = Lists.newArrayList();
        for (CopyNumberAlteration alteration : alterations) {
            final String cna = new StringJoiner(DELIMITER).add(circosContig(alteration.chromosome()))
                    .add(String.valueOf(alteration.start()))
                    .add(String.valueOf(alteration.end()))
                    .add(String.valueOf(alteration.minorAllelePloidy() - 1))
                    .toString();
            result.add(cna);
        }

        return result;
    }

    @NotNull
    private List<String> createScatter(@NotNull final List<Segment> segments, @NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (Segment segment : segments) {

            final String colorOption = ChainColor.color(segment.chainId());

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
            final boolean isStartFoldback = Links.findLink(startPosition, links).filter(Link::startFoldback).isPresent();
            String startGlyph = isStartFoldback ? "glyph=triangle,glyph_size=20" : "glyph=circle";
            if (segment.openStart()) {
                startGlyph = "glyph=square";
            }

            final StringJoiner start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                    .add(String.valueOf(segment.start()))
                    .add(String.valueOf(segment.start()))
                    .add(String.valueOf(segment.track()))
                    .add(colorOption + "," + startGlyph);
            result.add(start.toString());

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
            final boolean isEndFoldback = Links.findLink(endPosition, links).filter(Link::startFoldback).isPresent();
            String endGlyph = isEndFoldback ? "glyph=triangle,glyph_size=20" : "glyph=circle";
            if (segment.openEnd()) {
                endGlyph = "glyph=square";
            }

            final StringJoiner end = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                    .add(String.valueOf(segment.end()))
                    .add(String.valueOf(segment.end()))
                    .add(String.valueOf(segment.track()))
                    .add(colorOption + "," + endGlyph);
            result.add(end.toString());

        }

        return result;
    }

    @NotNull
    private List<String> createLinks(@NotNull final List<Segment> segments, @NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (final Link link : Links.clean(links)) {

            long linkUsage = segments.stream()
                    .filter(x -> x.chromosome().equals(link.startChromosome()) && (x.start() == link.startPosition()
                            || x.end() == link.startPosition()))
                    .count();

            final String linkString = new StringJoiner(DELIMITER).add(circosContig(link.startChromosome()))
                    .add(String.valueOf(link.startPosition()))
                    .add(String.valueOf(link.startPosition()))
                    .add(circosContig(link.endChromosome()))
                    .add(String.valueOf(link.endPosition()))
                    .add(String.valueOf(link.endPosition()))
                    .add(ChainColor.color(link.chainId()) + "," + thickness(linkUsage))
                    .toString();
            result.add(linkString);
        }

        return result;
    }

    @NotNull
    private List<String> createConnectors(int maxTracks, @NotNull final List<Segment> segments, @NotNull final List<Link> link) {
        final List<String> result = Lists.newArrayList();
        for (Segment segment : segments) {

            double r1 = CircosConfigWriter.svTrackPixels(maxTracks, segment.track());

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
            if (Links.findStartLink(startPosition, link).isPresent() || Links.findEndLink(startPosition, link).isPresent()) {
                long connectorUsage = segments.stream()
                        .filter(x -> x.chromosome().equals(segment.chromosome()) && x.start() == segment.start()
                                && x.track() >= segment.track())
                        .count();

                final String start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.start()))
                        .add("r1=" + r1 + "p," + ChainColor.color(segment.chainId()) + "," + thickness(connectorUsage))
                        .toString();
                result.add(start);
            }

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
            if (Links.findStartLink(endPosition, link).isPresent() || Links.findEndLink(endPosition, link).isPresent()) {
                long connectorUsage = segments.stream()
                        .filter(x -> x.chromosome().equals(segment.chromosome()) && x.end() == segment.end()
                                && x.track() >= segment.track())
                        .count();
                final String end = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.end()))
                        .add("r1=" + r1 + "p," + ChainColor.color(segment.chainId()) + "," + thickness(connectorUsage))
                        .toString();
                result.add(end);
            }

        }

        return result;
    }

    @NotNull
    private List<String> createKaryotypes(@NotNull final Map<String, Integer> contigLengths) {
        final List<String> result = Lists.newArrayList();
        for (String contig : contigLengths.keySet()) {

            final String start = new StringJoiner(" ").add("chr -")
                    .add(circosContig(contig))
                    .add(HumanChromosome.fromString(contig).toString())
                    .add(String.valueOf(1))
                    .add(String.valueOf(contigLengths.get(contig)))
                    .add("chr" + HumanChromosome.fromString(contig).toString())
                    .toString();
            result.add(start);
        }

        return result;
    }

    @NotNull
    private List<String> createHistogramTrack(@NotNull final Map<String, Integer> contigLengths, @NotNull final List<Segment> segments) {

        final List<String> result = Lists.newArrayList();
        for (Segment scaled : segments) {

            final int contigLength = contigLengths.get(scaled.chromosome());

            final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(1))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.track()))
                    .add("thickness=0")
                    .toString();
            result.add(start);

            final String entry = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.valueOf(scaled.track()))
                    .add(ChainColor.color(scaled.chainId()))
                    .toString();
            result.add(entry);

            final String end = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.valueOf(contigLength))
                    .add(String.valueOf(scaled.track()))
                    .add("thickness=0")
                    .toString();
            result.add(end);

        }

        return result;
    }

    @NotNull
    private Map<String, Integer> contigLengths(@NotNull final List<GenomePosition> positions) {
        final Map<String, Integer> results = Maps.newHashMap();
        for (GenomePosition position : positions) {
            int end = (int) position.position();
            results.merge(position.chromosome(), end, Math::max);
        }
        return results;
    }

    @NotNull
    private List<String> createPositionText(@NotNull final List<GenomePosition> originalLinks,
            @NotNull final List<GenomePosition> scaledLinks) {

        final Map<String, Integer> contigLengths = contigLengths(scaledLinks);
        final Set<String> result = Sets.newHashSet();

        for (int i = 0; i < originalLinks.size(); i++) {

            final GenomePosition original = originalLinks.get(i);
            final GenomePosition scaled = scaledLinks.get(i);
            if (scaled.position() != 1 && scaled.position() != contigLengths.get(scaled.chromosome())) {

                final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                        .add(String.valueOf(scaled.position()))
                        .add(String.valueOf(scaled.position()))
                        .add(String.format("%,d", original.position()))
                        .toString();

                result.add(start);
            }
        }

        return result.stream().sorted().distinct().collect(Collectors.toList());
    }

    @NotNull
    private static String circosContig(@NotNull final String chromosome) {
        return "hs" + HumanChromosome.fromString(chromosome);
    }

    @NotNull
    private static String thickness(long usage) {
        return "thickness=" + Math.max(1, (4 + (usage - 1) * 4));
    }

}
