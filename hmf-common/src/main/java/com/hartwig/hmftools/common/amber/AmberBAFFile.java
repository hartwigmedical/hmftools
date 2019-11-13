package com.hartwig.hmftools.common.amber;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class AmberBAFFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final Logger LOGGER = LogManager.getLogger(AmberBAFFile.class);

    private static final String DELIMITER = "\t";
    private static final String AMBER_EXTENSION = ".amber.baf.tsv";
    private static final String AMBER_EXTENSION_OLD = ".amber.baf";

    private AmberBAFFile() {
    }

    @NotNull
    public static String generateAmberFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    @NotNull
    public static String generateAmberFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + AMBER_EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + AMBER_EXTENSION_OLD;
    }

    @NotNull
    public static Multimap<Chromosome, AmberBAF> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull final Multimap<String, AmberBAF> bafs) throws IOException {
        List<AmberBAF> sortedBafs = Lists.newArrayList(bafs.values());
        Collections.sort(sortedBafs);
        write(filename, sortedBafs);
    }

    public static void write(@NotNull final String filename, @NotNull final List<AmberBAF> bafs) throws IOException {
        Files.write(new File(filename).toPath(), toLines(bafs));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<AmberBAF> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(AmberBAFFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome")
                .add("position")
                .add("tumorBAF")
                .add("tumorModifiedBAF")
                .add("tumorDepth")
                .add("normalBAF")
                .add("normalModifiedBAF")
                .add("normalDepth")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final AmberBAF ratio) {
        return new StringJoiner(DELIMITER).add(ratio.chromosome())
                .add(String.valueOf(ratio.position()))
                .add(FORMAT.format(ratio.tumorBAF()))
                .add(FORMAT.format(ratio.tumorModifiedBAF()))
                .add(String.valueOf(ratio.tumorDepth()))
                .add(Doubles.isFinite(ratio.normalBAF()) ? FORMAT.format(ratio.normalBAF()) : "0")
                .add(Doubles.isFinite(ratio.normalModifiedBAF()) ? FORMAT.format(ratio.normalModifiedBAF()) : "0")
                .add(String.valueOf(ratio.normalDepth()))
                .toString();
    }

    @NotNull
    private static Multimap<Chromosome, AmberBAF> fromLines(@NotNull final List<String> lines) {
        Multimap<Chromosome, AmberBAF> result = ArrayListMultimap.create();
        for (int i = 1; i < lines.size(); i++) {
            final String line = lines.get(i);
            try {
                final AmberBAF region = fromString(line);
                result.put(HumanChromosome.fromString(region.chromosome()), region);
            } catch (RuntimeException e) {
                LOGGER.info("Unable to parse line {}: {}", i, line);
                throw e;
            }
        }

        return result;
    }

    @NotNull
    private static AmberBAF fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        ImmutableAmberBAF.Builder builder = ImmutableAmberBAF.builder()
                .chromosome(values[0])
                .position(Long.parseLong(values[1]))
                .tumorBAF(Double.parseDouble(values[2]))
                .tumorDepth(0)
                .normalBAF(0.5)
                .normalDepth(0);

        if (values.length == 8) {
            builder.tumorDepth(Integer.parseInt(values[4]))
                    .normalBAF(Double.parseDouble(values[5]))
                    .normalDepth(Integer.parseInt(values[7]));
        }

        return builder.build();
    }
}
