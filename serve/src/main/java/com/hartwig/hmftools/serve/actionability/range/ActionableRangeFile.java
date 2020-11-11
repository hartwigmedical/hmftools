package com.hartwig.hmftools.serve.actionability.range;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEventFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActionableRangeFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_RANGE_TSV = "actionableRanges.tsv";

    private ActionableRangeFile() {
    }

    @NotNull
    public static String actionableRangeTsvPath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_RANGE_TSV;
    }

    public static void write(@NotNull String actionableRangeTsv, @NotNull List<ActionableRange> actionableRanges) throws IOException {
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
        //TODO: rangeInfo and feature can be removed when range positions are verified
        return new StringJoiner(DELIMITER).add("gene")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("range")
                .add("mutationType")
                .add("source")
                .add("treatment")
                .add("cancerType")
                .add("doid")
                .add("level")
                .add("direction")
                .add("url")
                .add("feature")
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
        String[] values = line.split(DELIMITER);
        String url = values.length > 11 ? values[11] : Strings.EMPTY;

        return ImmutableActionableRange.builder()
                .gene(values[0])
                .chromosome(values[1])
                .start(Long.parseLong(values[2]))
                .end(Long.parseLong(values[3]))
                .rangeInfo(0) //TODO: @Lieke, (changes to values[4] currently your to and from lines are inconsistent and files on datastore do not contain rangeInfo
                .mutationType(MutationTypeFilter.valueOf(values[4]))
                .source(ActionableEventFactory.sourceFromFileValue(values[5]))
                .treatment(values[6])
                .cancerType(values[7])
                .doid(values[8])
                .level(EvidenceLevel.valueOf(values[9]))
                .direction(ActionableEventFactory.directionFromFileValue(values[10]))
                .url(url)
                .feature(Lists.newArrayList()) //TODO: take feature list from file
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<ActionableRange> actionableRanges) {
        List<String> lines = Lists.newArrayList();
        for (ActionableRange actionableRange : actionableRanges) {
            lines.add(toLine(actionableRange));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull ActionableRange range) {
        //TODO: rangeInfo and feature can be removed when range positions are verified
        return new StringJoiner(DELIMITER).add(range.gene())
                .add(range.chromosome())
                .add(Long.toString(range.start()))
                .add(Long.toString(range.end()))
                .add(Integer.toString(range.rangeInfo()))
                .add(range.mutationType().toString())
                .add(range.source().display())
                .add(range.treatment())
                .add(range.cancerType())
                .add(range.doid())
                .add(range.level().toString())
                .add(range.direction().display())
                .add(range.url())
                .add(range.feature().toString())
                .toString();
    }
}
