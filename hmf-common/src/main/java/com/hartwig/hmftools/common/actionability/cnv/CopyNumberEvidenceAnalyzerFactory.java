package com.hartwig.hmftools.common.actionability.cnv;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class CopyNumberEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private CopyNumberEvidenceAnalyzerFactory() {
    }

    @NotNull
    public static CopyNumberEvidenceAnalyzer loadFromFileCNVs(String fileCNVs) throws IOException {
        final List<ActionableCopyNumber> CNVs = new ArrayList<>();
        final List<String> lineCNVs = Files.readAllLines(new File(fileCNVs).toPath());

        for (String line : lineCNVs) {
            if (!line.contains("event") || !line.contains("actionability")) {
                CNVs.add(fromLineCNVs(line));
            }
        }
        return new CopyNumberEvidenceAnalyzer(CNVs);
    }

    @NotNull
    private static ActionableCopyNumber fromLineCNVs(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableCopyNumber.builder()
                .gene(values[0])
                .cnvType(values[1])
                .source(values[2])
                .reference(values[3])
                .drug(values[4])
                .drugsType(values[5])
                .cancerType(values[6])
                .level(values[8])
                .response(values[11])
                .build();
    }
}
