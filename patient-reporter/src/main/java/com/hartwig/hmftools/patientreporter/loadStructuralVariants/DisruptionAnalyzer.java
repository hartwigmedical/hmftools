package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class DisruptionAnalyzer {

    @NotNull
    private final List<Disruption> disruption;

    DisruptionAnalyzer(@NotNull final List<Disruption> disruption) {
        this.disruption = disruption;
    }

    @NotNull
    public List<Disruption> filteringDisruptions() {
        List<Disruption> raportableDisruptions = Lists.newArrayList();
        for (Disruption disruptionItem: disruption) {
            if (disruptionItem.reportable()) {
                raportableDisruptions.add(disruptionItem);
            }
        }
        return raportableDisruptions;
    }
}
