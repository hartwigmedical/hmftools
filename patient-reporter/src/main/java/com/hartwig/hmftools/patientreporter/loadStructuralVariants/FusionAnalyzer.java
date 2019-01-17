package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionAnalyzer {

    @NotNull
    private final List<FusionReaderFile> fusionReaderFile;

    FusionAnalyzer(@NotNull final List<FusionReaderFile> fusionReaderFile) {
        this.fusionReaderFile = fusionReaderFile;
    }

    @Nullable
    public List<FusionReaderFile> filteringFusions() {
        List<FusionReaderFile> raportableFusions = Lists.newArrayList();
        for (int i = 1; i <= fusionReaderFile.size(); i ++) {
        }
        return raportableFusions;
    }
}
