package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;
import java.util.Map;

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
        Map<String, String> reportableGenesCanonicalTranscriptMap = geneModel.disruptionGeneCanonicalTranscriptMap();

        for (Disruption disruption: disruptions) {
            if (reportableGenesCanonicalTranscriptMap.keySet().contains(disruption.gene())) {
                if (reportableGenesCanonicalTranscriptMap.get(disruption.gene()).equals(disruption.transcript())) {
                    reportableDisruptions.add(disruption);
                }
            }
        }
        return reportableDisruptions;
    }
}
