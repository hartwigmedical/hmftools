package com.hartwig.hmftools.serve.actionability.range;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class RangeEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private RangeEvidenceAnalyzerFactory() {

    }

    @NotNull
    public static RangeEvidenceAnalyzer loadFromFileRange(@NotNull String actionableRangeTsv) throws IOException {
        final List<ActionableRange> ranges = Lists.newArrayList();
        final List<String> lineRange = Files.readAllLines(new File(actionableRangeTsv).toPath());

        // Skip header line for range
        for (String lineRanges : lineRange.subList(1, lineRange.size())) {
            ranges.add(fromLineRange(lineRanges));
        }
        return new RangeEvidenceAnalyzer(ranges);
    }

    @NotNull
    private static ActionableRange fromLineRange(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableRange.builder()
                .gene(values[0])
                .chromosome(values[1])
                .start(values[2])
                .end(values[3])
                .mutationType(values[4])
                .source(values[5])
                .drug(values[6])
                .drugType(values[7])
                .cancerType(values[8])
                .level(values[9])
                .direction(values[10])
                .link(values[11])
                .build();
    }
}
