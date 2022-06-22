package com.hartwig.hmftools.serve.extraction.copynumber;

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
import com.hartwig.hmftools.serve.extraction.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;

public final class KnownCopyNumberFile {

    private static final String DELIMITER = "\t";
    private static final String KNOWN_COPY_NUMBER_TSV = "KnownCopyNumbers.SERVE.tsv";

    private KnownCopyNumberFile() {
    }

    @NotNull
    public static String knownCopyNumberTsvPath(@NotNull String outputDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(outputDir + File.separator + KNOWN_COPY_NUMBER_TSV);
    }

    @NotNull
    public static List<KnownCopyNumber> read(@NotNull String file) throws IOException {
        List<String> lines = Files.readAllLines(new File(file).toPath());

        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    @VisibleForTesting
    static List<KnownCopyNumber> fromLines(@NotNull List<String> lines) {
        List<KnownCopyNumber> copyNumber = Lists.newArrayList();
        for (String line : lines) {
            copyNumber.add(fromLine(line));
        }
        return copyNumber;
    }

    @NotNull
    private static KnownCopyNumber fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableKnownCopyNumber.builder()
                .gene(values[0])
                .type(CopyNumberType.valueOf(values[1]))
                .sources(Knowledgebase.fromCommaSeparatedSourceString(values[2]))
                .build();
    }

    public static void write(@NotNull String copyNumberTsv, @NotNull Iterable<KnownCopyNumber> copyNumbers) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(copyNumbers));
        Files.write(new File(copyNumberTsv).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("gene").add("type").add("sources").toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull Iterable<KnownCopyNumber> copyNumbers) {
        List<String> lines = Lists.newArrayList();
        for (KnownCopyNumber copyNumber : sort(copyNumbers)) {
            lines.add(toLine(copyNumber));
        }
        return lines;
    }

    @NotNull
    private static List<KnownCopyNumber> sort(@NotNull Iterable<KnownCopyNumber> copyNumbers) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<KnownCopyNumber> sorted = Lists.newArrayList(copyNumbers);
        sorted.sort(new KnownCopyNumberComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull KnownCopyNumber copyNumber) {
        return new StringJoiner(DELIMITER).add(copyNumber.gene())
                .add(copyNumber.type().toString().toLowerCase())
                .add(Knowledgebase.toCommaSeparatedSourceString(copyNumber.sources()))
                .toString();
    }
}
