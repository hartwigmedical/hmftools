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
    public static FusionEvidenceAnalyzer loadFromFileFusion(@NotNull String actionableFusionTsv) throws IOException {
        final List<ActionableFusion> fusions = Lists.newArrayList();
        final List<String> lineFusion = Files.readAllLines(new File(actionableFusionTsv).toPath());

        // Skip header line for fusions
        for (String lineFusions : lineFusion.subList(1, lineFusion.size())) {
            fusions.add(fromLineFusion(lineFusions));
        }
        return new FusionEvidenceAnalyzer(fusions);
    }

    @NotNull
    private static ActionableFusion fromLineFusion(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableFusion.builder()
                .fusion(values[0])
                .source(values[1])
                .drug(values[2])
                .drugType(values[3])
                .cancerType(values[4])
                .level(values[5])
                .direction(values[6])
                .link(values[7])
                .build();
    }
}
