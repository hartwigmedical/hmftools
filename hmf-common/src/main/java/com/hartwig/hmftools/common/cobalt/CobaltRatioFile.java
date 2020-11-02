package com.hartwig.hmftools.common.cobalt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public final class CobaltRatioFile {
    public static final String TUMOR_ONLY_REFERENCE_SAMPLE = "DIPLOID";

    private static final DecimalFormat FORMAT = new DecimalFormat("#.####");

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".cobalt.ratio.tsv";
    private static final String EXTENSION_OLD = ".cobalt";

    private CobaltRatioFile() {
    }

    @NotNull
    public static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + EXTENSION_OLD;
    }

    @NotNull
    public static ListMultimap<Chromosome, CobaltRatio> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull Multimap<Chromosome, CobaltRatio> ratios) throws IOException {
        List<CobaltRatio> sorted = Lists.newArrayList(ratios.values());
        Collections.sort(sorted);
        write(fileName, sorted);
    }

    private static void write(@NotNull final String fileName, @NotNull List<CobaltRatio> ratios) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(ratios));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<CobaltRatio> ratio) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(CobaltRatioFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome")
                .add("position")
                .add("referenceReadCount")
                .add("tumorReadCount")
                .add("referenceGCRatio")
                .add("tumorGCRatio")
                .add("referenceGCDiploidRatio")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final CobaltRatio position) {
        return new StringJoiner(DELIMITER).add(position.chromosome())
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.referenceReadCount()))
                .add(String.valueOf(position.tumorReadCount()))
                .add(String.valueOf(FORMAT.format(position.referenceGCRatio())))
                .add(String.valueOf(FORMAT.format(position.tumorGCRatio())))
                .add(String.valueOf(FORMAT.format(position.referenceGCDiploidRatio())))
                .toString();
    }

    @NotNull
    private static ListMultimap<Chromosome, CobaltRatio> fromLines(@NotNull final List<String> lines) {
        final ListMultimap<Chromosome, CobaltRatio> result = ArrayListMultimap.create();
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i);
            final CobaltRatio ratio = fromLine(line);
            result.put(HumanChromosome.fromString(ratio.chromosome()), ratio);
        }
        return result;
    }

    @NotNull
    private static CobaltRatio fromLine(@NotNull final String ratioLine) {
        final String[] values = ratioLine.split(DELIMITER);

        final String chromosome = values[0].trim();
        final long position = Long.parseLong(values[1].trim());

        return ImmutableCobaltRatio.builder()
                .chromosome(chromosome)
                .position(position)
                .referenceReadCount(Integer.parseInt(values[2].trim()))
                .tumorReadCount(Integer.parseInt(values[3].trim()))
                .referenceGCRatio(Double.parseDouble(values[4].trim()))
                .tumorGCRatio(Double.parseDouble(values[5].trim()))
                .referenceGCDiploidRatio(Double.parseDouble(values[6].trim()))
                .build();
    }
}
