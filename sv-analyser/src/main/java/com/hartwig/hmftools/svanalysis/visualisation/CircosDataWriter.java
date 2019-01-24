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
        Files.write(new File(textPath).toPath(), createPositionText(unadjustedSegments, segments));

        final String histogramPath = filePrefix + ".histogram.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(contigLengths, segments));

        final String karyotypePath = filePrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(contigLengths));

        final String connectorPath = filePrefix + ".connector.circos";
        Files.write(new File(connectorPath).toPath(), createConnectors(maxTracks, segments, links));

        final String linkPath = filePrefix + ".link.circos";
        Files.write(new File(linkPath).toPath(), createLinks(links));

        final String scatterPath = filePrefix + ".scatter.circos";
        Files.write(new File(scatterPath).toPath(), createScatter(segments, links));

        final String cnaPath = filePrefix + ".cna.circos";
        Files.write(new File(cnaPath).toPath(), createCNA(alterations));

        final String mapPath = filePrefix + ".map.circos";
        Files.write(new File(mapPath).toPath(), createMinorAllelePloidy(alterations));

        final String terminals = filePrefix + ".terminals.circos";
        Files.write(new File(terminals).toPath(), createTerminals(segments));

    }

    @NotNull
    private List<String> createTerminals(@NotNull final List<Segment> segments) {
        final List<String> result = Lists.newArrayList();
        for (final Segment segment : segments) {
            if (segment.startTerminal() != SegmentTerminal.NONE) {
                final String cna = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.track() + (segment.startTerminal() == SegmentTerminal.TELOMERE ? "T" : "C")))
                        .add(ChainColor.color(segment.chainId()))
                        .toString();
                result.add(cna);
            }

            if (segment.endTerminal() != SegmentTerminal.NONE) {
                final String cna = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.track() + (segment.endTerminal() == SegmentTerminal.TELOMERE ? "T" : "C")))
                        .add(ChainColor.color(segment.chainId()))
                        .toString();
                result.add(cna);
            }
        }

        return result;
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
            final boolean isStartFoldback =
                    Links.findLink(startPosition, links).filter(x -> x.startType() == Link.Type.FOLDBACK).isPresent();
            String startGlyph = isStartFoldback ? "glyph=triangle,glyph_size=20" : "glyph=circle";
            if (segment.startTerminal() != SegmentTerminal.NONE) {
                startGlyph = "glyph=square,glyph_size=1";
            }

            final StringJoiner start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                    .add(String.valueOf(segment.start()))
                    .add(String.valueOf(segment.start()))
                    .add(String.valueOf(segment.track()))
                    .add(colorOption + "," + startGlyph);
            result.add(start.toString());

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
            final boolean isEndFoldback = Links.findLink(endPosition, links).filter(x -> x.endType() == Link.Type.FOLDBACK).isPresent();
            String endGlyph = isEndFoldback ? "glyph=triangle,glyph_size=20" : "glyph=circle";
            if (segment.endTerminal() != SegmentTerminal.NONE) {
                endGlyph = "glyph=square,glyph_size=1";
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
    private List<String> createLinks(@NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (final Link link : Links.clean(links)) {

            final String linkString = new StringJoiner(DELIMITER).add(circosContig(link.startChromosome()))
                    .add(String.valueOf(link.startPosition()))
                    .add(String.valueOf(link.startPosition()))
                    .add(circosContig(link.endChromosome()))
                    .add(String.valueOf(link.endPosition()))
                    .add(String.valueOf(link.endPosition()))
                    .add(ChainColor.color(link.chainId()) + "," + thickness(link.traverseCount()))
                    .toString();
            result.add(linkString);
        }

        return result;
    }

    @NotNull
    private List<String> createConnectors(int maxTracks, @NotNull final List<Segment> segments, @NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (Segment segment : segments) {

            double r1 = CircosConfigWriter.svTrackPixels(maxTracks, segment.track());

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
            int startLinkUsage = Links.linkTraverseCount(startPosition, links);

            if (startLinkUsage > 0) {
                long segmentsBelow = segments.stream()
                        .filter(x -> x.chromosome().equals(segment.chromosome()) && x.start() == segment.start()
                                && x.track() < segment.track())
                        .count();

                final String start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.start()))
                        .add("r1=" + r1 + "p," + ChainColor.color(segment.chainId()) + "," + thickness(startLinkUsage - segmentsBelow))
                        .toString();
                result.add(start);
            }

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
            int endLinkUsage = Links.linkTraverseCount(endPosition, links);

            if (endLinkUsage > 0) {
                long segmentsBelow = segments.stream()
                        .filter(x -> x.chromosome().equals(segment.chromosome()) && x.end() == segment.end() && x.track() < segment.track())
                        .count();
                final String end = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.end()))
                        .add("r1=" + r1 + "p," + ChainColor.color(segment.chainId()) + "," + thickness(endLinkUsage - segmentsBelow))
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
                    .add(ChainColor.color(scaled.chainId()) + "," + thickness(scaled.traverseCount()))
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
    private List<String> createPositionText(@NotNull final List<Segment> originalLinks, @NotNull final List<Segment> scaledLinks) {

        final Set<String> result = Sets.newHashSet();

        for (int i = 0; i < originalLinks.size(); i++) {

            final Segment original = originalLinks.get(i);
            final Segment scaled = scaledLinks.get(i);

            String startText = original.startTerminal() == SegmentTerminal.CENTROMERE ? "Centromere" : "Telomere";
            startText = original.startTerminal() == SegmentTerminal.NONE ? String.format("%,d", original.start()) : startText;

            final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.start()))
                    .add(startText)
                    .toString();

            result.add(start);

            String endText = original.endTerminal() == SegmentTerminal.CENTROMERE ? "Centromere" : "Telomere";
            endText = original.endTerminal() == SegmentTerminal.NONE ? String.format("%,d", original.start()) : endText;

            final String end = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.valueOf(scaled.end()))
                    .add(endText)
                    .toString();

            result.add(end);
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
