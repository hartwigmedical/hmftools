package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class FusionAnalyzer {

    @NotNull
    private final List<FusionReaderFile> fusionReaderFile;

    FusionAnalyzer(@NotNull final List<FusionReaderFile> fusionReaderFile) {
        this.fusionReaderFile = fusionReaderFile;
    }

    @NotNull
    public List<FusionReaderFile> filteringFusions() {
        List<FusionReaderFile> raportableFusions = Lists.newArrayList();
        for (FusionReaderFile fusion: fusionReaderFile){
            if (fusion.reportable().equals(true)) {
                raportableFusions.add(fusion);
            }
        }
        return raportableFusions;
    }
}
