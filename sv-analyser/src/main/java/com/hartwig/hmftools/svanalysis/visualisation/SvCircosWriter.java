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

public class SvCircosWriter {

    private static final String BACKGROUND_COLOUR = "246,251,244";
    private static final String DELIMITER = "\t";
    private static final int POSITION_BUFFER = 2;

    private final String outputPrefix;

    public SvCircosWriter(@NotNull final String outputPrefix) {
        this.outputPrefix = outputPrefix;
    }

    public void writeLinks(@NotNull final List<GenomeRegion> rawRegions, @NotNull final List<SvLink> rawLinks) throws IOException {

        final ScalePosition scaler = new ScalePosition(1 + POSITION_BUFFER, rawRegions);

        final List<SvRegion> scaledRegions = scaledLinks(scaler.scaleRegions(rawRegions));
        final List<SvLink> scaledLinks = scaler.scaleLinks(rawLinks);

        final String textPath = outputPrefix + ".text.circos";
        Files.write(new File(textPath).toPath(), createPositionText(rawRegions, scaledRegions));

        final String histogramPath = outputPrefix + ".histogram.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(scaledRegions));

        final String karyotypePath = outputPrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(scaledRegions));

        final String connectorPath = outputPrefix + ".connector.circos";
        Files.write(new File(connectorPath).toPath(), createConnectors(0.3, 0.6/10, scaledRegions));

        final String linkPath = outputPrefix + ".link.circos";
        Files.write(new File(linkPath).toPath(), createLinks(scaledLinks));

    }

    @NotNull
    private List<String> createLinks(@NotNull final List<SvLink> svLinks) {
        final List<String> result = Lists.newArrayList();
        for (SvLink svLink : svLinks) {
            final String link = new StringJoiner(DELIMITER)
                    .add(circosContig(svLink.startChromosome()))
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
    private List<String> createConnectors(double r0, double radiusChange, @NotNull final List<SvRegion> scaledLinks) {
        final List<String> result = Lists.newArrayList();
        for (SvRegion scaled : scaledLinks) {

            double r1 = r0 + (scaled.value() - 1) * radiusChange;

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
    private List<String> createKaryotypes(@NotNull final List<SvRegion> scaledLinks) {
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
    private List<String> createHistogramTrack(@NotNull final List<SvRegion> scaledLinks) {
        final Map<String, Integer> contigLengths = contigLengths(scaledLinks);

        final List<String> result = Lists.newArrayList();
        for (SvRegion scaled : scaledLinks) {

            final int contigLength = contigLengths.get(scaled.chromosome());

            final String start = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(1))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.value()))
                    .add("color=" + BACKGROUND_COLOUR)
                    .toString();
            result.add(start);

            final String entry = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.start()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.valueOf(scaled.value()))
                    .toString();
            result.add(entry);

            final String end = new StringJoiner(DELIMITER).add(circosContig(scaled.chromosome()))
                    .add(String.valueOf(scaled.end()))
                    .add(String.valueOf(contigLength))
                    .add(String.valueOf(scaled.value()))
                    .add("color=" + BACKGROUND_COLOUR)
                    .toString();
            result.add(end);

        }

        return result;
    }

    @NotNull
    private Map<String, Integer> contigLengths(@NotNull final List<SvRegion> scaledLinks) {
        final Map<String, Integer> results = Maps.newHashMap();
        for (SvRegion scaledLink : scaledLinks) {
            int end = (int) scaledLink.end() + POSITION_BUFFER;
            results.merge(scaledLink.chromosome(), end, Math::max);
        }
        return results;
    }

    @NotNull
    private List<String> createPositionText(@NotNull final List<GenomeRegion> originalLinks, @NotNull final List<SvRegion> scaledLinks) {

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

    @NotNull
    private List<SvRegion> scaledLinks(@NotNull final List<GenomeRegion> scaledLinks) {

        final Map<String, Integer> valueMap = Maps.newHashMap();
        final List<SvRegion> result = Lists.newArrayList();

        int currentValue = 1;
        for (GenomeRegion scaledLink : scaledLinks) {

            if (!valueMap.containsKey(scaledLink.chromosome())) {
                valueMap.put(scaledLink.chromosome(), currentValue);
            } else {
                currentValue = Math.max(currentValue, valueMap.get(scaledLink.chromosome()) + 1);
                valueMap.put(scaledLink.chromosome(), currentValue);
            }

            result.add(new SvRegion(scaledLink, currentValue));
        }

        return result;

    }

}
