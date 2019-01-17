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
        for (int i = 1; i < fusionReaderFile.size(); i ++) {
            if (fusionReaderFile.get(i).reportable().equals(true)){
                raportableFusions.add(fusionReaderFile.get(i));
            }
        }
        return raportableFusions;
    }
}
