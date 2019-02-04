package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class FusionAnalyzer {

    @NotNull
    private final List<Fusion> fusions;

    FusionAnalyzer(@NotNull final List<Fusion> fusions) {
        this.fusions = fusions;
    }

    @NotNull
    @VisibleForTesting
    List<Fusion> fusions() {
        return fusions;
    }

    @NotNull
    public List<Fusion> reportableFusions() {
        List<Fusion> reportableFusions = Lists.newArrayList();
        for (Fusion fusion: fusions){
            if (fusion.reportable()) {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }
}
