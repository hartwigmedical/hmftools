package com.hartwig.hmftools.serve.extraction.exon;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;
import com.hartwig.hmftools.serve.util.RefGenomeVersion;

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

    public static List<KnownExon> read(@NotNull String file) throws IOException {
        List<String> lines = Files.readAllLines(new File(file).toPath());

        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static List<KnownExon> fromLines(@NotNull List<String> lines) {
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
                        .chromosome(values[1])
                        .start(Long.parseLong(values[2]))
                        .end(Long.parseLong(values[3]))
                        .mutationType(MutationTypeFilter.valueOf(values[4]))
                        .exonEnsemblId(values[5])
                        .exonIndex(Integer.parseInt(values[6]))
                        .build())
                .sources(Sets.newHashSet(Knowledgebase.DOCM)) //TODO fix
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
                .add("chromosome")
                .add("start")
                .add("end")
                .add("mutationType")
                .add("exonEnsemblId")
                .add("exonIndex")
                .add("sources")
                .toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull Iterable<KnownExon> exons) {
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
                .add(exon.annotation().chromosome())
                .add(String.valueOf(exon.annotation().start()))
                .add(String.valueOf(exon.annotation().end()))
                .add(exon.annotation().mutationType().toString())
                .add(exon.annotation().exonEnsemblId())
                .add(String.valueOf(exon.annotation().exonIndex()))
                .add(Knowledgebase.commaSeparatedSourceString(exon.sources()))
                .toString();
    }
}
