package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class DisruptionAnalyzer {

    @NotNull
    private final List<DisruptionReaderFile> disruptionReaderFiles;

    DisruptionAnalyzer(@NotNull final List<DisruptionReaderFile> disruptionReaderFiles) {
        this.disruptionReaderFiles = disruptionReaderFiles;
    }

}
