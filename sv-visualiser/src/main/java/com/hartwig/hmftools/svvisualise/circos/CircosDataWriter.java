package com.hartwig.hmftools.svvisualise.circos;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.LinkedHashMap;
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
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.svvisualise.data.CopyNumberAlteration;
import com.hartwig.hmftools.svvisualise.data.Exon;
import com.hartwig.hmftools.svvisualise.data.Link;
import com.hartwig.hmftools.svvisualise.data.Links;
import com.hartwig.hmftools.svvisualise.data.Segment;

import org.jetbrains.annotations.NotNull;

public class CircosDataWriter {

    private static final int MIN_KAROTYPE_LENGTH = 10;
    private static final String DELIMITER = "\t";

    private final ColorPicker colorPicker;
    private final String filePrefix;
    private final int maxTracks;
    private boolean debug;

    public CircosDataWriter(final boolean debug, final ColorPicker colorPicker, @NotNull final String sample,
            @NotNull final String outputDir, final int maxTracks) {
        this.debug = debug;
        this.colorPicker = colorPicker;
        this.filePrefix = outputDir + File.separator + sample;
        this.maxTracks = maxTracks;
    }

    public void write(@NotNull final List<Segment> unadjustedSegments, @NotNull final List<Link> unadjustedLinks,
            @NotNull final List<CopyNumberAlteration> unadjustedAlterations, @NotNull final List<Exon> unadjustedExons) throws IOException {

        final List<GenomeRegion> unadjustedFragileSites =
                Highlights.limitHighlightsToSegments(Highlights.fragileSites(), unadjustedSegments);

        final List<GenomeRegion> unadjustedLineElements =
                Highlights.limitHighlightsToSegments(Highlights.lineElements(), unadjustedSegments);

        // Note we do not add exons here because we want them interpolated.
        final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
        unadjustedPositions.addAll(Links.allPositions(unadjustedLinks));
        unadjustedPositions.addAll(Span.allPositions(unadjustedSegments));
        unadjustedPositions.addAll(Span.allPositions(unadjustedAlterations));
        unadjustedPositions.addAll(Span.allPositions(unadjustedFragileSites));
        unadjustedPositions.addAll(Span.allPositions(unadjustedLineElements));

        final ScalePosition scalePosition = new ScalePosition(unadjustedPositions);
        final List<GenomePosition> scaledPositions = scalePosition.scaled();
        final Map<String, Integer> contigLengths = contigLengths(scaledPositions);

        final List<Segment> segments = scalePosition.scaleSegments(unadjustedSegments);
        final List<Link> links = scalePosition.scaleLinks(unadjustedLinks);
        final List<CopyNumberAlteration> alterations = scalePosition.scaleAlterations(unadjustedAlterations);
        final List<GenomeRegion> fragileSites = scalePosition.scaleRegions(unadjustedFragileSites);
        final List<GenomeRegion> lineElements = scalePosition.scaleRegions(unadjustedLineElements);
        final List<Exon> exons = scalePosition.interpolateExons(unadjustedExons);

        final String exonPath = filePrefix + ".exon.circos";
        Files.write(new File(exonPath).toPath(), exons(exons));

        final String genePath = filePrefix + ".gene.circos";
        Files.write(new File(genePath).toPath(), genes(exons));

        final String geneNamePath = filePrefix + ".gene.name.circos";
        Files.write(new File(geneNamePath).toPath(), geneName(exons));

        final String textPath = filePrefix + ".text.circos";
        Files.write(new File(textPath).toPath(), createPositionText(debug, unadjustedLinks, links, segments));

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

        final String fragile = filePrefix + ".fragile.circos";
        Files.write(new File(fragile).toPath(), highlights(fragileSites));

        final String line = filePrefix + ".line_element.circos";
        Files.write(new File(line).toPath(), highlights(lineElements));

        final String distances = filePrefix + ".distance.circos";
        if (unadjustedAlterations.size() < 200) {
            Files.write(new File(distances).toPath(), createDistances(unadjustedAlterations, alterations));
        } else {
            Files.write(new File(distances).toPath(), Collections.emptySet());
        }

    }

    @NotNull
    private List<String> genes(@NotNull final List<Exon> exons) {
        final List<String> result = Lists.newArrayList();
        final Set<String> genes = exons.stream().map(Exon::gene).collect(Collectors.toSet());
        for (final String gene : genes) {
            final List<Exon> geneExons = exons.stream().filter(x -> x.gene().equals(gene)).collect(toList());
            long min = geneExons.stream().mapToLong(GenomeRegion::start).min().orElse(0);
            long max = geneExons.stream().mapToLong(GenomeRegion::end).max().orElse(0);

            final String exonString = new StringJoiner(DELIMITER).add(circosContig(geneExons.get(0).chromosome()))
                    .add(String.valueOf(min))
                    .add(String.valueOf(max))
                    .add(String.valueOf(1))
                    .toString();
            result.add(exonString);

        }

        return result;
    }

    @NotNull
    private List<String> geneName(@NotNull final List<Exon> exons) {
        final List<String> result = Lists.newArrayList();
        final Set<String> genes = exons.stream().map(Exon::gene).collect(Collectors.toSet());
        for (final String gene : genes) {
            final List<Exon> geneExons = exons.stream().filter(x -> x.gene().equals(gene)).collect(toList());
            long min = geneExons.stream().mapToLong(GenomeRegion::start).min().orElse(0);

            final String geneName = geneExons.get(0).gene();
            final double labelSize = geneNameLabelSize(geneName);

            final String exonString = new StringJoiner(DELIMITER).add(circosContig(geneExons.get(0).chromosome()))
                    .add(String.valueOf(min))
                    .add(String.valueOf(min))
                    .add(geneName)
                    .add("label_size=" +labelSize + "p,rpadding=0r")
                    .toString();
            result.add(exonString);
        }

        return result;
    }

    private static double geneNameLabelSize(@NotNull final String gene) {
        double availablePixels = CircosConfigWriter.PIXELS * (CircosConfigWriter.EXON_OUTER_RADIUS - CircosConfigWriter.EXON_INNER_RADIUS);
        return Math.min(26, 4 + Math.floor(availablePixels / gene.length()));
    }

    @NotNull
    private List<String> exons(@NotNull final List<Exon> exons) {
        final List<String> result = Lists.newArrayList();
        final Set<String> chromosome = exons.stream().map(GenomeRegion::chromosome).collect(Collectors.toSet());
        for (String contig : chromosome) {

            final List<GenomeRegion> contigRegions = exons.stream()
                    .filter(x -> x.chromosome().equals(contig))
                    .map(x -> GenomeRegionFactory.create(x.chromosome(), x.start(), x.end()))
                    .sorted()
                    .distinct()
                    .collect(Collectors.toList());

            for (int i = 0; i < contigRegions.size(); i++) {
                final GenomeRegion region = contigRegions.get(i);
                final String exonString = new StringJoiner(DELIMITER).add(circosContig(region.chromosome()))
                        .add(String.valueOf(region.start()))
                        .add(String.valueOf(region.end()))
                        .add(String.valueOf(1))
                        .toString();
                result.add(exonString);

                if (i < contigRegions.size() - 1) {
                    final GenomeRegion next = contigRegions.get(i + 1);

                    final String betweenString = new StringJoiner(DELIMITER).add(circosContig(region.chromosome()))
                            .add(String.valueOf(region.end()))
                            .add(String.valueOf(next.start()))
                            .add(String.valueOf(0))
                            .toString();
                    //                    result.add(betweenString);
                }
            }
        }

        return result;
    }

    @NotNull
    private List<String> highlights(@NotNull final List<GenomeRegion> regions) {
        return regions.stream()
                .map(x -> new StringJoiner(DELIMITER).add(circosContig(x.chromosome()))
                        .add(String.valueOf(x.start()))
                        .add(String.valueOf(x.end()))
                        .toString())
                .collect(toList());
    }

    @NotNull
    private List<String> createDistances(@NotNull final List<CopyNumberAlteration> unadjustedSegment,
            @NotNull final List<CopyNumberAlteration> segments) {

        final Map<String, Integer> contigLengths = contigLengthsFromRegions(segments);

        final List<String> result = Lists.newArrayList();
        for (int i = 0; i < unadjustedSegment.size(); i++) {
            final CopyNumberAlteration adjusted = segments.get(i);
            final CopyNumberAlteration unadjusted = unadjustedSegment.get(i);

            if (adjusted.start() != 1 && adjusted.end() != contigLengths.get(adjusted.chromosome())) {
                final String distance = new StringJoiner(DELIMITER).add(circosContig(adjusted.chromosome()))
                        .add(String.valueOf(adjusted.start()))
                        .add(String.valueOf(adjusted.end()))
                        .add(shorthand(unadjusted.end() - unadjusted.start()))
                        .toString();
                result.add(distance);
            }
        }
        return result;
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
                        .add(colorPicker.color(segment.clusterId(), segment.chainId()))
                        .toString();
                result.add(cna);
            }

            if (segment.endTerminal() != SegmentTerminal.NONE) {
                final String cna = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.track() + (segment.endTerminal() == SegmentTerminal.TELOMERE ? "T" : "C")))
                        .add(colorPicker.color(segment.clusterId(), segment.chainId()))
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

            final String colorOption = colorPicker.color(segment.clusterId(), segment.chainId());

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
            final boolean isStartFoldback =
                    Links.findStartLink(startPosition, links).filter(x -> x.startInfo().equals("FOLDBACK")).isPresent()
                            || Links.findEndLink(startPosition, links).filter(x -> x.endInfo().equals("FOLDBACK")).isPresent();
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
            final boolean isEndFoldback =
                    Links.findStartLink(endPosition, links).filter(x -> x.startInfo().equals("FOLDBACK")).isPresent() || Links.findEndLink(
                            endPosition,
                            links).filter(x -> x.endInfo().equals("FOLDBACK")).isPresent();

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
                    .add(colorPicker.color(link.clusterId(), link.chainId()) + "," + thickness(link.traverseCount()))
                    .toString();
            result.add(linkString);
        }

        return result;
    }

    @NotNull
    private List<String> createConnectors(int maxTracks, @NotNull final List<Segment> segments, @NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (Segment segment : segments) {

            final double r1 = CircosConfigWriter.svTrackPixels(maxTracks, segment.track());

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
                        .add("r1=" + r1 + "p," + colorPicker.connectorColor(segment.clusterId(), segment.chainId()) + "," + thickness(
                                startLinkUsage - segmentsBelow))
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
                        .add("r1=" + r1 + "p," + colorPicker.connectorColor(segment.clusterId(), segment.chainId()) + "," + thickness(
                                endLinkUsage - segmentsBelow))
                        .toString();
                result.add(end);
            }

        }

        double rTrack1 = CircosConfigWriter.svTrackPixels(maxTracks, 0);
        for (Link link : links) {
            if (link.isSimpleSV()) {
                if (link.isValidStart()) {
                    final String start = new StringJoiner(DELIMITER).add(circosContig(link.startChromosome()))
                            .add(String.valueOf(link.startPosition()))
                            .add(String.valueOf(link.startPosition()))
                            .add("r1=" + rTrack1 + "p," + colorPicker.connectorColor(link.clusterId(), link.chainId()) + "," + thickness(
                                    link.traverseCount()))
                            .toString();
                    result.add(start);
                }

                if (link.isValidEnd()) {
                    final String end = new StringJoiner(DELIMITER).add(circosContig(link.endChromosome()))
                            .add(String.valueOf(link.endPosition()))
                            .add(String.valueOf(link.endPosition()))
                            .add("r1=" + rTrack1 + "p," + colorPicker.connectorColor(link.clusterId(), link.chainId()) + "," + thickness(
                                    link.traverseCount()))
                            .toString();
                    result.add(end);
                }

            }

        }

        return result.stream().sorted().collect(toList());
    }

    @NotNull
    private List<String> createKaryotypes(@NotNull final Map<String, Integer> contigLengths) {
        final List<String> result = Lists.newArrayList();
        for (String contig : contigLengths.keySet()) {

            final String start = new StringJoiner(" ").add("chr -")
                    .add(circosContig(contig))
                    .add(HumanChromosome.fromString(contig).toString())
                    .add(String.valueOf(1))
                    .add(String.valueOf(Math.max(MIN_KAROTYPE_LENGTH, contigLengths.get(contig))))
                    .add("chr" + HumanChromosome.fromString(contig).toString())
                    .toString();
            result.add(start);
        }

        return result;
    }

    @NotNull
    private List<String> createHistogramTrack(@NotNull final Map<String, Integer> contigLengths, @NotNull final List<Segment> segments) {

        final List<String> result = Lists.newArrayList();
        for (Segment segment : segments) {

            final int contigLength = contigLengths.get(segment.chromosome());

            final String start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                    .add(String.valueOf(1))
                    .add(String.valueOf(segment.start()))
                    .add(String.valueOf(segment.track()))
                    .add("thickness=0")
                    .toString();
            result.add(start);

            final String entry = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                    .add(String.valueOf(segment.start()))
                    .add(String.valueOf(segment.end()))
                    .add(String.valueOf(segment.track()))
                    .add(colorPicker.color(segment.clusterId(), segment.chainId()) + "," + thickness(segment.traverseCount()))
                    .toString();
            result.add(entry);

            final String end = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                    .add(String.valueOf(segment.end()))
                    .add(String.valueOf(contigLength))
                    .add(String.valueOf(segment.track()))
                    .add("thickness=0")
                    .toString();
            result.add(end);

        }

        return result;
    }

    @NotNull
    private Map<String, Integer> contigLengths(@NotNull final List<GenomePosition> positions) {
        final Map<String, Integer> results = new LinkedHashMap<>();
        final List<GenomePosition> sortedPositions = positions.stream().sorted().collect(toList());

        for (GenomePosition position : sortedPositions) {
            int end = (int) position.position();
            results.merge(position.chromosome(), end, Math::max);
        }
        return results;
    }

    @NotNull
    private Map<String, Integer> contigLengthsFromRegions(@NotNull final List<? extends GenomeRegion> positions) {
        final Map<String, Integer> results = Maps.newHashMap();
        for (GenomeRegion region : positions) {
            int end = (int) region.end();
            results.merge(region.chromosome(), end, Math::max);
        }
        return results;
    }

    @NotNull
    private List<String> createPositionText(boolean debug, @NotNull final List<Link> originalLinks, @NotNull final List<Link> scaledLinks,
            @NotNull final List<Segment> segments) {

        final Set<String> result = Sets.newHashSet();

        for (int i = 0; i < originalLinks.size(); i++) {

            final Link original = originalLinks.get(i);
            final Link scaled = scaledLinks.get(i);

            if (scaled.isValidStart()) {
                final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.startChromosome()))
                        .add(String.valueOf(scaled.startPosition()))
                        .add(String.valueOf(scaled.startPosition()))
                        .add(String.valueOf(debug ? original.svId() : original.startPosition()))
                        .toString();

                result.add(start);
            }

            if (scaled.isValidEnd()) {
                final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.endChromosome()))
                        .add(String.valueOf(scaled.endPosition()))
                        .add(String.valueOf(scaled.endPosition()))
                        .add(String.valueOf(debug ? original.svId() : original.endPosition()))
                        .toString();

                result.add(start);
            }
        }

        for (final Segment segment : segments) {
            if (segment.startTerminal() != SegmentTerminal.NONE) {
                final String startText = segment.startTerminal() == SegmentTerminal.CENTROMERE ? "Centromere" : "Telomere";
                final String start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.start()))
                        .add(startText)
                        .toString();
                result.add(start);
            }

            if (segment.endTerminal() != SegmentTerminal.NONE) {
                final String endText = segment.endTerminal() == SegmentTerminal.CENTROMERE ? "Centromere" : "Telomere";
                final String start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.end()))
                        .add(endText)
                        .toString();
                result.add(start);
            }

        }

        return result.stream().sorted().distinct().collect(toList());
    }

    @NotNull
    private static String circosContig(@NotNull final String chromosome) {
        return "hs" + HumanChromosome.fromString(chromosome);
    }

    @NotNull
    private static String thickness(long usage) {
        return "thickness=" + Math.max(1, Math.min(10, usage) * 4);
    }

    @NotNull
    static String shorthand(long value) {
        if (value < 100) {
            return String.valueOf(value);
        }

        if (value < 99_950) {
            return String.format("%.1fk", value / 1_000d);
        }

        return String.format("%.1fm", value / 1_000_000d);
    }

}
