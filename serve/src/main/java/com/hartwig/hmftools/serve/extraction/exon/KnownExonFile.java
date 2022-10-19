package com.hartwig.hmftools.serve.extraction.exon;

import static com.hartwig.hmftools.common.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.datamodel.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;

public final class KnownExonFile {

    private static final String DELIMITER = "\t";
    private static final String KNOWN_EXON_TSV = "KnownExons.SERVE.tsv";

    private KnownExonFile() {
    }

    @NotNull
    public static String knownExonTsvPath(@NotNull String outputDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(outputDir + File.separator + KNOWN_EXON_TSV);
    }

    @NotNull
    public static List<KnownExon> read(@NotNull String file) throws IOException {
        List<String> lines = Files.readAllLines(new File(file).toPath());

        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    @VisibleForTesting
    static List<KnownExon> fromLines(@NotNull List<String> lines) {
        List<KnownExon> exons = Lists.newArrayList();
        for (String line : lines) {
            exons.add(fromLine(line));
        }
        return exons;
    }

    @NotNull
    private static KnownExon fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableKnownExon.builder()
                .annotation(ImmutableExonAnnotation.builder()
                        .gene(values[0])
                        .transcript(values[1])
                        .chromosome(values[2])
                        .start(Integer.parseInt(values[3]))
                        .end(Integer.parseInt(values[4]))
                        .mutationType(MutationTypeFilter.valueOf(values[5]))
                        .rank(Integer.parseInt(values[6]))
                        .build())
                .sources(Knowledgebase.fromCommaSeparatedSourceString(values[7]))
                .build();
    }

    public static void write(@NotNull String exonTsv, @NotNull Iterable<KnownExon> exons) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(exons));
        Files.write(new File(exonTsv).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("gene")
                .add("transcript")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("mutationType")
                .add("exonRank")
                .add("sources")
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<KnownExon> exons) {
        List<String> lines = Lists.newArrayList();
        for (KnownExon exon : sort(exons)) {
            lines.add(toLine(exon));
        }
        return lines;
    }

    @NotNull
    private static List<KnownExon> sort(@NotNull Iterable<KnownExon> codons) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<KnownExon> sorted = Lists.newArrayList(codons);
        sorted.sort(new KnownExonComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull KnownExon exon) {
        return new StringJoiner(DELIMITER).add(exon.annotation().gene())
                .add(exon.annotation().transcript())
                .add(exon.annotation().chromosome())
                .add(String.valueOf(exon.annotation().start()))
                .add(String.valueOf(exon.annotation().end()))
                .add(exon.annotation().mutationType().toString())
                .add(String.valueOf(exon.annotation().rank()))
                .add(Knowledgebase.toCommaSeparatedSourceString(exon.sources()))
                .toString();
    }
}
