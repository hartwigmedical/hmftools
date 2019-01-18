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

    public void writeLinks(@NotNull final List<GenomeRegion> rawLinks) throws IOException {

        final List<SvLink> scaledLinks = scaledLinks(rawLinks);

        final String textPath = outputPrefix + ".text.circos";
        Files.write(new File(textPath).toPath(), createPositionText(rawLinks, scaledLinks));

        final String histogramPath = outputPrefix + ".histogram.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(scaledLinks));

        final String karyotypePath = outputPrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(scaledLinks));

    }

    private List<String> createKaryotypes(@NotNull final List<SvLink> scaledLinks) {
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

    private List<String> createHistogramTrack(@NotNull final List<SvLink> scaledLinks) {
        final Map<String, Integer> contigLengths = contigLengths(scaledLinks);

        final List<String> result = Lists.newArrayList();
        for (SvLink scaled : scaledLinks) {

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

    private Map<String, Integer> contigLengths(@NotNull final List<SvLink> scaledLinks) {
        final Map<String, Integer> results = Maps.newHashMap();
        for (SvLink scaledLink : scaledLinks) {
            int end = (int) scaledLink.end() + POSITION_BUFFER;
            results.merge(scaledLink.chromosome(), end, Math::max);
        }
        return results;
    }

    private List<String> createPositionText(@NotNull final List<GenomeRegion> originalLinks, @NotNull final List<SvLink> scaledLinks) {

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

    private List<SvLink> scaledLinks(@NotNull final List<GenomeRegion> rawLinks) {

        final Map<String, Integer> valueMap = Maps.newHashMap();
        final List<SvLink> result = Lists.newArrayList();
        final List<GenomeRegion> scaledLinks = ScalePosition.scale(1 + POSITION_BUFFER, rawLinks);

        int currentValue = 1;
        for (GenomeRegion scaledLink : scaledLinks) {

            if (!valueMap.containsKey(scaledLink.chromosome())) {
                valueMap.put(scaledLink.chromosome(), currentValue);
            } else {
                currentValue = Math.max(currentValue, valueMap.get(scaledLink.chromosome()) + 1);
                valueMap.put(scaledLink.chromosome(), currentValue);
            }

            result.add(new SvLink(scaledLink, currentValue));
        }

        return result;

    }

}
