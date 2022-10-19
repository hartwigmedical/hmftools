package com.hartwig.hmftools.serve.extraction.fusion;

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
import com.hartwig.hmftools.common.serve.datamodel.fusion.FusionPairComparator;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class KnownFusionPairFile {

    private static final String DELIMITER = "\t";
    private static final String KNOWN_FUSION_PAIR_TSV = "KnownFusionPairs.SERVE.tsv";

    private KnownFusionPairFile() {
    }

    @NotNull
    public static String knownFusionPairTsvPath(@NotNull String outputDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(outputDir + File.separator + KNOWN_FUSION_PAIR_TSV);
    }

    @NotNull
    public static List<KnownFusionPair> read(@NotNull String file) throws IOException {
        List<String> lines = Files.readAllLines(new File(file).toPath());

        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    @VisibleForTesting
    static List<KnownFusionPair> fromLines(@NotNull List<String> lines) {
        List<KnownFusionPair> fusionPairs = Lists.newArrayList();
        for (String line : lines) {
            fusionPairs.add(fromLine(line));
        }
        return fusionPairs;
    }

    @NotNull
    private static KnownFusionPair fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableKnownFusionPair.builder()
                .geneUp(values[0])
                .minExonUp(values[1].isEmpty() ? null: Integer.valueOf(values[1]))
                .maxExonUp(values[2].isEmpty() ? null: Integer.valueOf(values[2]))
                .geneDown(values[3])
                .minExonDown(values[4].isEmpty() ? null: Integer.valueOf(values[4]))
                .maxExonDown(values[5].isEmpty() ? null: Integer.valueOf(values[5]))
                .sources(Knowledgebase.fromCommaSeparatedSourceString(values[6]))
                .build();
    }

    public static void write(@NotNull String fusionPairTsv, @NotNull Iterable<KnownFusionPair> fusionPairs) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(fusionPairs));
        Files.write(new File(fusionPairTsv).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("geneUp")
                .add("minExonUp")
                .add("maxExonUp")
                .add("geneDown")
                .add("minExonDown")
                .add("maxExonDown")
                .add("sources")
                .toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull Iterable<KnownFusionPair> fusionPairs) {
        List<String> lines = Lists.newArrayList();
        for (KnownFusionPair fusionPair : sort(fusionPairs)) {
            lines.add(toLine(fusionPair));
        }
        return lines;
    }

    @NotNull
    private static List<KnownFusionPair> sort(@NotNull Iterable<KnownFusionPair> fusionPairs) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<KnownFusionPair> sorted = Lists.newArrayList(fusionPairs);
        sorted.sort(new FusionPairComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull KnownFusionPair fusionPair) {
        return new StringJoiner(DELIMITER).add(fusionPair.geneUp())
                .add(nullToEmpty(fusionPair.minExonUp()))
                .add(nullToEmpty(fusionPair.maxExonUp()))
                .add(fusionPair.geneDown())
                .add(nullToEmpty(fusionPair.minExonDown()))
                .add(nullToEmpty(fusionPair.maxExonDown()))
                .add(Knowledgebase.toCommaSeparatedSourceString(fusionPair.sources()))
                .toString();
    }

    @NotNull
    private static String nullToEmpty(@Nullable Integer integer) {
        return integer != null ? String.valueOf(integer) : Strings.EMPTY;
    }
}
