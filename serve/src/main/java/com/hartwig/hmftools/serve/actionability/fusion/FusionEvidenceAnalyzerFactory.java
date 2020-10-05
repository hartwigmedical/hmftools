package com.hartwig.hmftools.serve.actionability.fusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class FusionEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private FusionEvidenceAnalyzerFactory() {
    }

    @NotNull
    public static FusionEvidenceAnalyzer loadFromActionableFusionTsv(@NotNull String actionableFusionTsv) throws IOException {
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(actionableFusionTsv).toPath());

        // Skip header line for fusions
        for (String line : lines.subList(1, lines.size())) {
            actionableFusions.add(fromLine(line));
        }
        return new FusionEvidenceAnalyzer(actionableFusions);
    }

    @NotNull
    private static ActionableFusion fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableActionableFusion.builder()
                .fusion(values[0])
                .source(values[1])
                .drug(values[2])
                .cancerType(values[3])
                .doid(values[4])
                .level(values[5])
                .direction(values[6])
                .build();
    }
}
