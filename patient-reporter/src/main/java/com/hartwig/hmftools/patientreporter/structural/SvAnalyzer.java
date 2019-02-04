package com.hartwig.hmftools.patientreporter.structural;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SvAnalyzer {

    @NotNull
    private final List<Fusion> fusions;
    @NotNull
    private final List<Disruption> disruptions;

    @NotNull
    public static SvAnalyzer fromFiles(@NotNull String fusionFile, @NotNull String disruptionFile) throws IOException {
        List<Fusion> fusions = FusionFileReader.fromFusionFile(fusionFile);
        List<Disruption> disruptions = DisruptionFileReader.fromDisruptionFile(disruptionFile);
        return new SvAnalyzer(fusions, disruptions);
    }

    private SvAnalyzer(@NotNull final List<Fusion> fusions, @NotNull final List<Disruption> disruptions) {
        this.fusions = fusions;
        this.disruptions = disruptions;
    }

    @NotNull
    public SvAnalysis run(@NotNull GeneModel geneModel, @NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        List<ReportableGeneFusion> reportableFusions = ReportableGeneFusionFactory.convert(reportableFusions());
        List<ReportableGeneDisruption> reportableGeneDisruptions =
                ReportableGeneDisruptionFactory.convert(reportableDisruptions(geneModel), geneCopyNumbers);

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<SimpleGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(toSimpleGeneFusions(reportableFusions), primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerFusion);

        return ImmutableSvAnalysis.builder()
                .reportableFusions(reportableFusions)
                .reportableDisruptions(reportableGeneDisruptions)
                .evidenceItems(filteredEvidence)
                .build();
    }

    @NotNull
    private List<Fusion> reportableFusions() {
        List<Fusion> reportableFusions = Lists.newArrayList();
        for (Fusion fusion : fusions) {
            if (fusion.reportable()) {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }

    @NotNull
    private List<Disruption> reportableDisruptions(@NotNull GeneModel geneModel) {
        List<Disruption> reportableDisruptions = Lists.newArrayList();
        Set<String> reportableGenes = geneModel.disruptionGenePanel();

        for (Disruption disruption : disruptions) {
            if (reportableGenes.contains(disruption.gene()) && disruption.canonical() && disruption.isDisruptive()) {
                reportableDisruptions.add(disruption);
            }
        }

        return reportableDisruptions;
    }

    @NotNull
    private static List<SimpleGeneFusion> toSimpleGeneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        List<SimpleGeneFusion> simpleGeneFusions = Lists.newArrayList();
        for (ReportableGeneFusion fusionReport : fusions) {
            simpleGeneFusions.add(ImmutableSimpleGeneFusion.builder()
                    .fiveGene(fusionReport.geneStart())
                    .threeGene(fusionReport.geneEnd())
                    .build());
        }
        return simpleGeneFusions;
    }
}
