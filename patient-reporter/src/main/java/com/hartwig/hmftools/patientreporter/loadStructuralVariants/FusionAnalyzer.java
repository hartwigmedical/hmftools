package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class FusionAnalyzer {

    @NotNull
    private final List<Fusion> fusion;

    FusionAnalyzer(@NotNull final List<Fusion> fusion) {
        this.fusion = fusion;
    }

    @NotNull
    public List<Fusion> filteringFusions() {
        List<Fusion> reportableFusions = Lists.newArrayList();
        for (Fusion fusion: fusion){
            if (fusion.reportable()) {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }
}
