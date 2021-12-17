package com.hartwig.hmftools.serve.extraction.codon;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;

public final class KnownCodonFile {

    private static final String DELIMITER = "\t";
    private static final String KNOWN_CODON_TSV = "KnownCodons.SERVE.tsv";

    private KnownCodonFile() {
    }

    @NotNull
    public static String knownCodonTsvPath(@NotNull String outputDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(outputDir + File.separator + KNOWN_CODON_TSV);
    }

    @NotNull
    public static List<KnownCodon> read(@NotNull String file) throws IOException {
        List<String> lines = Files.readAllLines(new File(file).toPath());

        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    @VisibleForTesting
    static List<KnownCodon> fromLines(@NotNull List<String> lines) {
        List<KnownCodon> codons = Lists.newArrayList();
        for (String line : lines) {
            codons.add(fromLine(line));
        }
        return codons;
    }

    @NotNull
    private static KnownCodon fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableKnownCodon.builder()
                .annotation(ImmutableCodonAnnotation.builder()
                        .gene(values[0])
                        .transcript(values[1])
                        .chromosome(values[2])
                        .start(Integer.parseInt(values[3]))
                        .end(Integer.parseInt(values[4]))
                        .mutationType(MutationTypeFilter.valueOf(values[5]))
                        .rank(Integer.parseInt(values[6]))
                        .build())
                .sources(Knowledgebase.fromCommaSeparatedTechnicalDisplayString(values[7]))
                .build();
    }

    public static void write(@NotNull String codonTsv, @NotNull Iterable<KnownCodon> codons) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(codons));
        Files.write(new File(codonTsv).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("gene")
                .add("transcript")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("mutationType")
                .add("codonRank")
                .add("sources")
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<KnownCodon> codons) {
        List<String> lines = Lists.newArrayList();
        for (KnownCodon codon : sort(codons)) {
            lines.add(toLine(codon));
        }
        return lines;
    }

    @NotNull
    private static List<KnownCodon> sort(@NotNull Iterable<KnownCodon> codons) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<KnownCodon> sorted = Lists.newArrayList(codons);
        sorted.sort(new KnownCodonComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull KnownCodon codon) {
        return new StringJoiner(DELIMITER).add(codon.annotation().gene())
                .add(codon.annotation().transcript())
                .add(codon.annotation().chromosome())
                .add(String.valueOf(codon.annotation().start()))
                .add(String.valueOf(codon.annotation().end()))
                .add(codon.annotation().mutationType().toString())
                .add(String.valueOf(codon.annotation().rank()))
                .add(Knowledgebase.toCommaSeparatedSourceString(codon.sources()))
                .toString();
    }
}
