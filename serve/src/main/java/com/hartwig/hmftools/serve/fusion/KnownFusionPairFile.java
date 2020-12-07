package com.hartwig.hmftools.serve.fusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class KnownFusionPairFile {

    private static final String DELIMITER = "\t";

    private KnownFusionPairFile() {
    }

    public static void write(@NotNull String fusionPairTsv, @NotNull List<KnownFusionPair> fusionPairs) throws IOException {
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
                .add("geneUp")
                .add("sources")
                .toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull List<KnownFusionPair> fusionPairs) {
        List<String> lines = Lists.newArrayList();
        for (KnownFusionPair fusionPair : sort(fusionPairs)) {
            lines.add(toLine(fusionPair));
        }
        return lines;
    }

    @NotNull
    @VisibleForTesting
    static List<KnownFusionPair> sort(@NotNull List<KnownFusionPair> fusionPairs) {
        // Need to make a copy since the input list may be immutable and cannot be sorted!
        List<KnownFusionPair> sorted = Lists.newArrayList(fusionPairs);
        sorted.sort(new KnownFusionPairComparator());

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
                .add(Knowledgebase.commaSeparatedSourceString(fusionPair.sources()))
                .toString();
    }

    @NotNull
    private static String nullToEmpty(@Nullable Integer integer) {
        return integer != null ? String.valueOf(integer) : Strings.EMPTY;
    }
}
