package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;

public class DisruptionAnalyzer {

    @NotNull
    private final List<Disruption> disruptions;

    DisruptionAnalyzer(@NotNull final List<Disruption> disruptions) {
        this.disruptions = disruptions;
    }

    @NotNull
    public List<Disruption> reportableDisruptions(@NotNull GeneModel geneModel) {
        List<Disruption> reportableDisruptions = Lists.newArrayList();
        Set<String> reportableGenes = geneModel.disruptionGenePanel();

        for (Disruption disruption : disruptions) {
            if (reportableGenes.contains(disruption.gene()) && disruption.canonical() && disruption.isDisruptive()) {
                reportableDisruptions.add(disruption);
            }
        }

        return reportableDisruptions;
    }
}
