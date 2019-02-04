package com.hartwig.hmftools.patientreporter.structural;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;

public class SvAnalyzer {

    @NotNull
    private final List<Fusion> fusions;
    @NotNull
    private final List<Disruption> disruptions;

    @NotNull
    public static SvAnalyzer fromFiles(@NotNull String fusionFile, @NotNull String disruptionFile) throws IOException {
        List<Fusion> fusions = FusionFactory.fromFusionFile(fusionFile);
        List<Disruption> disruptions = DisruptionFactory.fromDisruptionFile(disruptionFile);
        return new SvAnalyzer(fusions, disruptions);
    }

    private SvAnalyzer(@NotNull final List<Fusion> fusions, @NotNull final List<Disruption> disruptions) {
        this.fusions = fusions;
        this.disruptions = disruptions;
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
