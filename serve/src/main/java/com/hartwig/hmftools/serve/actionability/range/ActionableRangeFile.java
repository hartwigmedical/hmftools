package com.hartwig.hmftools.serve.actionability.range;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;

public final class ActionableRangeFile {

    private static final String ACTIONABLE_RANGE_TSV = "ActionableRanges.tsv";

    private ActionableRangeFile() {
    }

    @NotNull
    public static String actionableRangeTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_RANGE_TSV);
    }

    public static void write(@NotNull String actionableRangeTsv, @NotNull Iterable<ActionableRange> actionableRanges) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableRanges));
        Files.write(new File(actionableRangeTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableRange> read(@NotNull String actionableRangeTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableRangeTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("gene")
                .add("transcript")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("mutationType")
                .add("rangeType")
                .add("rank")
                .add(ActionableFileFunctions.header())
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<ActionableRange> fromLines(@NotNull List<String> lines) {
        List<ActionableRange> actionableRanges = Lists.newArrayList();
        for (String line : lines) {
            actionableRanges.add(fromLine(line));
        }
        return actionableRanges;
    }

    @NotNull
    private static ActionableRange fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActionableRange.builder()
                .from(ActionableFileFunctions.fromLine(values, 8))
                .gene(values[0])
                .transcript(values[1])
                .chromosome(values[2])
                .start(Integer.parseInt(values[3]))
                .end(Integer.parseInt(values[4]))
                .mutationType(MutationTypeFilter.valueOf(values[5]))
                .rangeType(RangeType.valueOf(values[6]))
                .rank(Integer.parseInt(values[7]))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<ActionableRange> actionableRanges) {
        List<String> lines = Lists.newArrayList();
        for (ActionableRange actionableRange : sort(actionableRanges)) {
            lines.add(toLine(actionableRange));
        }
        return lines;
    }

    @NotNull
    private static List<ActionableRange> sort(@NotNull Iterable<ActionableRange> actionableRanges) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<ActionableRange> sorted = Lists.newArrayList(actionableRanges);
        sorted.sort(new ActionableRangeComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull ActionableRange range) {
        return new StringJoiner(FIELD_DELIMITER).add(range.gene())
                .add(range.transcript())
                .add(range.chromosome())
                .add(Long.toString(range.start()))
                .add(Long.toString(range.end()))
                .add(range.mutationType().toString())
                .add(range.rangeType().toString())
                .add(String.valueOf(range.rank()))
                .add(ActionableFileFunctions.toLine(range))
                .toString();
    }
}
