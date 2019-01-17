package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class FusionAnalyzer {

    @NotNull
    private final List<FusionReaderFile> fusionReaderFile;

    FusionAnalyzer(@NotNull final List<FusionReaderFile> fusionReaderFile) {
        this.fusionReaderFile = fusionReaderFile;
    }
}
