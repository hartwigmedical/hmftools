package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Optional;
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

import org.jetbrains.annotations.NotNull;

public class CircosDataWriter {

    private static final String BACKGROUND_COLOUR = "246,251,244";
    private static final String DELIMITER = "\t";
    private static final int POSITION_BUFFER = 2;

    private final String filePrefix;
    private final int maxTracks;

    public CircosDataWriter(@NotNull final String sample, @NotNull final String outputDir, final int maxTracks) {
        this.filePrefix = outputDir + File.separator + sample;
        this.maxTracks = maxTracks;
    }

    public void write(@NotNull final List<Track> unadjustedTracks, @NotNull final List<Link> unadjustedLinks,
            @NotNull final List<CopyNumberAlteration> unadjustedAlterations) throws IOException {

        final List<GenomeRegion> unadjustedRegions = Lists.newArrayList();
        unadjustedRegions.addAll(unadjustedTracks);
        unadjustedRegions.addAll(unadjustedAlterations);

        final ScalePosition scalePosition = new ScalePosition(POSITION_BUFFER, unadjustedRegions);
        final List<Track> tracks = scalePosition.scaleTracks(unadjustedTracks);
        final List<Link> links = scalePosition.scaleLinks(unadjustedLinks);
        final List<CopyNumberAlteration> alterations = scalePosition.scaleAlterations(unadjustedAlterations);

        final String textPath = filePrefix + ".text.circos";
        Files.write(new File(textPath).toPath(), createPositionText(scalePosition.original(), scalePosition.scaled()));

        final String histogramPath = filePrefix + ".histogram.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(tracks));

        final String karyotypePath = filePrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(tracks));

        final String connectorPath = filePrefix + ".connector.circos";
        Files.write(new File(connectorPath).toPath(), createConnectors(390, (1350d - 390d) / maxTracks, tracks, links));

        final String linkPath = filePrefix + ".link.circos";
        Files.write(new File(linkPath).toPath(), createLinks(links));

        final String scatterPath = filePrefix + ".scatter.circos";
        Files.write(new File(scatterPath).toPath(), createScatter(tracks, links));

        final String cnaPath = filePrefix + ".cna.circos";
        Files.write(new File(cnaPath).toPath(), createCNA(alterations));
    }

    @NotNull
    private List<String> createCNA(@NotNull final List<CopyNumberAlteration> alterations) {
        final List<String> result = Lists.newArrayList();
        for (CopyNumberAlteration alteration : alterations) {
            final String cna = new StringJoiner(DELIMITER).add(circosContig(alteration.chromosome()))
                    .add(String.valueOf(alteration.start()))
                    .add(String.valueOf(alteration.end()))
                    .add(String.valueOf(alteration.copyNumber()))
                    .toString();
            result.add(cna);
        }

        return result;
    }

    @NotNull
    private List<String> createScatter(@NotNull final List<Track> tracks, @NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (Track track : tracks) {

            final GenomePosition startPosition = GenomePositions.create(track.chromosome(), track.start());
            final boolean isStartFoldback =
                    startLink(startPosition, links).filter(Link::startFoldback).isPresent() || endLink(startPosition,
                            links).filter(Link::endFoldback).isPresent();

            final StringJoiner start = new StringJoiner(DELIMITER).add(circosContig(track.chromosome()))
                    .add(String.valueOf(track.start()))
                    .add(String.valueOf(track.start()))
                    .add(String.valueOf(track.track()))
                    .add(ChainColor.color(track.chainId()));
            //            if (isStartFoldback) {
            //                start.add("glyph=triangle");
            //            }

            result.add(start.toString());

            final GenomePosition endPosition = GenomePositions.create(track.chromosome(), track.end());
            final boolean isEndFoldback = startLink(endPosition, links).filter(Link::startFoldback).isPresent() || endLink(endPosition,
                    links).filter(Link::endFoldback).isPresent();
            final StringJoiner end = new StringJoiner(DELIMITER).add(circosContig(track.chromosome()))
                    .add(String.valueOf(track.end()))
                    .add(String.valueOf(track.end()))
                    .add(String.valueOf(track.track()))
                    .add(ChainColor.color(track.chainId()));
            //            if (isEndFoldback) {
            //                start.add("glyph=triangle");
            //            }
            result.add(end.toString());

        }

        return result;
    }

    @NotNull
    private List<String> createLinks(@NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (Link link : links) {
            final String linkString = new StringJoiner(DELIMITER).add(circosContig(link.startChromosome()))
                    .add(String.valueOf(link.startPosition()))
                    .add(String.valueOf(link.startPosition()))
                    .add(circosContig(link.endChromosome()))
                    .add(String.valueOf(link.endPosition()))
                    .add(String.valueOf(link.endPosition()))
                    .add(ChainColor.color(link.chainId()))
                    .toString();
            result.add(linkString);
        }

        return result;
    }

    @NotNull
    private List<String> createConnectors(double r0, double radiusChange, @NotNull final List<Track> tracks,
            @NotNull final List<Link> link) {
        final List<String> result = Lists.newArrayList();
        for (Track track : tracks) {

            double r1 = r0 + track.track() * radiusChange;

            final GenomePosition startPosition = GenomePositions.create(track.chromosome(), track.start());
            if (startLink(startPosition, link).isPresent() || endLink(startPosition, link).isPresent()) {
                final String start = new StringJoiner(DELIMITER).add(circosContig(track.chromosome()))
                        .add(String.valueOf(track.start()))
                        .add(String.valueOf(track.start()))
                        .add("r1=" + r1 + "p," + ChainColor.color(track.chainId()))
                        .toString();
                result.add(start);
            }

            final GenomePosition endPosition = GenomePositions.create(track.chromosome(), track.end());
            if (startLink(endPosition, link).isPresent() || endLink(endPosition, link).isPresent()) {
                final String end = new StringJoiner(DELIMITER).add(circosContig(track.chromosome()))
                        .add(String.valueOf(track.end()))
                        .add(String.valueOf(track.end()))
                        .add("r1=" + r1 + "p," + ChainColor.color(track.chainId()))
                        .toString();
                result.add(end);
            }

        }

        return result;
    }

    @NotNull
    private List<String> createKaryotypes(@NotNull final List<Track> tracks) {
        final Map<String, Integer> contigLengths = contigLengths(tracks);
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
    private List<String> createHistogramTrack(@NotNull final List<Track> tracks) {
        final Map<String, Integer> contigLengths = contigLengths(tracks);

        final List<String> result = Lists.newArrayList();
        for (Track scaled : tracks) {

            final int contigLength = contigLengths.get(scaled.chromosome());

            final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(1))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.track()))
                    .add("color=" + BACKGROUND_COLOUR)
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
                    .add("color=" + BACKGROUND_COLOUR)
                    .toString();
            result.add(end);

        }

        return result;
    }

    @NotNull
    private Map<String, Integer> contigLengths(@NotNull final List<Track> tracks) {
        final Map<String, Integer> results = Maps.newHashMap();
        for (Track scaledLink : tracks) {
            int end = (int) scaledLink.end() + POSITION_BUFFER;
            results.merge(scaledLink.chromosome(), end, Math::max);
        }
        return results;
    }

    @NotNull
    private List<String> createPositionText(@NotNull final List<GenomePosition> originalLinks,
            @NotNull final List<GenomePosition> scaledLinks) {

        Set<String> result = Sets.newHashSet();

        for (int i = 0; i < originalLinks.size(); i++) {
            final GenomePosition original = originalLinks.get(i);
            final GenomePosition scaled = scaledLinks.get(i);

            final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.position()))
                    .add(String.valueOf(scaled.position()))
                    .add(String.format("%,d", original.position()))
                    .toString();

            result.add(start);
        }

        return result.stream().sorted().distinct().collect(Collectors.toList());
    }

    @NotNull
    private static String circosContig(@NotNull final String chromosome) {
        return "hs" + HumanChromosome.fromString(chromosome);
    }

    static Optional<Link> endLink(@NotNull final GenomePosition position, @NotNull List<Link> links) {
        return links.stream()
                .filter(x -> x.endChromosome().equals(position.chromosome()) && x.endPosition() == position.position())
                .findFirst();
    }

    static Optional<Link> startLink(@NotNull final GenomePosition position, @NotNull List<Link> links) {
        return links.stream()
                .filter(x -> x.startChromosome().equals(position.chromosome()) && x.startPosition() == position.position())
                .findFirst();
    }

}
