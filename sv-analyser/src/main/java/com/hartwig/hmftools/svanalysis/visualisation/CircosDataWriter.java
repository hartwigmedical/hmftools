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

    public void write(@NotNull final List<Track> unadjustedTracks, @NotNull final List<Link> unadjustedLinks) throws IOException {

        final ScalePosition scalePosition = new ScalePosition(POSITION_BUFFER, unadjustedTracks);
        final List<Track> tracks = scalePosition.scaleTracks(unadjustedTracks);
        final List<Link> links = scalePosition.scaleLinks(unadjustedLinks);

        final String textPath = filePrefix + ".text.circos";
        Files.write(new File(textPath).toPath(), createPositionText(unadjustedTracks, tracks));

        final String histogramPath = filePrefix + ".histogram.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(tracks));

        final String karyotypePath = filePrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(tracks));

        final String connectorPath = filePrefix + ".connector.circos";
        Files.write(new File(connectorPath).toPath(), createConnectors(0.3, 0.6 / (maxTracks), tracks));

        final String linkPath = filePrefix + ".link.circos";
        Files.write(new File(linkPath).toPath(), createLinks(links));

    }

    @NotNull
    private List<String> createLinks(@NotNull final List<Link> links) {
        final List<String> result = Lists.newArrayList();
        for (Link svLink : links) {
            final String link = new StringJoiner(DELIMITER).add(circosContig(svLink.startChromosome()))
                    .add(String.valueOf(svLink.startPosition()))
                    .add(String.valueOf(svLink.startPosition()))
                    .add(circosContig(svLink.endChromosome()))
                    .add(String.valueOf(svLink.endPosition()))
                    .add(String.valueOf(svLink.endPosition()))
                    .toString();
            result.add(link);
        }

        return result;
    }

    @NotNull
    private List<String> createConnectors(double r0, double radiusChange, @NotNull final List<Track> scaledLinks) {
        final List<String> result = Lists.newArrayList();
        for (Track scaled : scaledLinks) {

            double r1 = r0 + (scaled.track() - 1) * radiusChange;

            final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.start()))
                    .add("r1=" + r1 + "r")
                    .toString();
            result.add(start);

            final String end = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.valueOf(scaled.end()))
                    .add("r1=" + r1 + "r")
                    .toString();
            result.add(end);

        }

        return result;
    }

    @NotNull
    private List<String> createKaryotypes(@NotNull final List<Track> scaledLinks) {
        final Map<String, Integer> contigLengths = contigLengths(scaledLinks);
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
    private List<String> createHistogramTrack(@NotNull final List<Track> scaledLinks) {
        final Map<String, Integer> contigLengths = contigLengths(scaledLinks);

        final List<String> result = Lists.newArrayList();
        for (Track scaled : scaledLinks) {

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
    private Map<String, Integer> contigLengths(@NotNull final List<Track> scaledLinks) {
        final Map<String, Integer> results = Maps.newHashMap();
        for (Track scaledLink : scaledLinks) {
            int end = (int) scaledLink.end() + POSITION_BUFFER;
            results.merge(scaledLink.chromosome(), end, Math::max);
        }
        return results;
    }

    @NotNull
    private List<String> createPositionText(@NotNull final List<Track> originalLinks, @NotNull final List<Track> scaledLinks) {

        Set<String> result = Sets.newHashSet();

        for (int i = 0; i < originalLinks.size(); i++) {
            final GenomeRegion original = originalLinks.get(i);
            final GenomeRegion scaled = scaledLinks.get(i);

            final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.start()))
                    .add(String.format("%,d", original.start()))
                    .toString();

            final String end = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.format("%,d", original.end()))
                    .toString();

            result.add(start);
            result.add(end);

        }

        return result.stream().sorted().distinct().collect(Collectors.toList());
    }

    @NotNull
    static String circosContig(@NotNull final String chromosome) {
        return "hs" + HumanChromosome.fromString(chromosome);
    }
}
