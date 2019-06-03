package com.hartwig.hmftools.patientreporter.structural;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SvAnalyzer {

    @NotNull
    private final List<ReportableGeneFusion> fusions;
    @NotNull
    private final List<ReportableDisruption> disruptions;

    @NotNull
    private final List<Fusion> fusionsOld;
    @NotNull
    private final List<Disruption> disruptionsOld;

    @NotNull
    public static SvAnalyzer fromFiles(@NotNull String fusionFile, @NotNull String disruptionFile) throws IOException
    {

        List<ReportableGeneFusion> fusions = Lists.newArrayList();
        List<ReportableDisruption> disruptions = Lists.newArrayList();

        List<Fusion> fusionsOld = Lists.newArrayList();
        List<Disruption> disruptionsOld = Lists.newArrayList();

        if(fusionFile.endsWith(ReportableGeneFusionFile.FILE_EXTENSION) && disruptionFile.endsWith(ReportableDisruptionFile.FILE_EXTENSION))
        {
            fusions = ReportableGeneFusionFile.read(fusionFile);
            disruptions = ReportableDisruptionFile.read(disruptionFile);

        }
        else
        {
            fusionsOld = FusionFileReader.fromFusionFile(fusionFile);
            disruptionsOld = DisruptionFileReader.fromDisruptionFile(disruptionFile);
        }

        return new SvAnalyzer(fusions, fusionsOld, disruptions, disruptionsOld);
    }

    @VisibleForTesting
    SvAnalyzer(@NotNull final List<ReportableGeneFusion> fusions, @NotNull final List<Fusion> fusionsOld,
            @NotNull final List<ReportableDisruption> disruptions, @NotNull final List<Disruption> disruptionsOld) {
        this.fusions = fusions;
        this.fusionsOld = fusionsOld;
        this.disruptions = disruptions;
        this.disruptionsOld = disruptionsOld;
    }

    @NotNull
    public SvAnalysis run(@NotNull GeneModel geneModel, @NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {

        List<ReportableGeneFusion> reportableFusions = !this.fusions.isEmpty()
                ? this.fusions : ReportableGeneFusionFactory.convert(reportableFusions());

        List<ReportableDisruption> reportableDisruptions = getReportableDisruptions(geneModel);

        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(reportableDisruptions, geneCopyNumbers);

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
        for (Fusion fusion : fusionsOld) {
            if (fusion.reportable()) {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }


    @NotNull
    private List<ReportableDisruption> getReportableDisruptions(@NotNull GeneModel geneModel)
    {
        Set<String> reportableGenes = geneModel.disruptionGenes();

        if(!disruptions.isEmpty())
        {
            return disruptions.stream()
                    .filter(x -> geneModel.disruptionGenes().contains(x.gene()))
                    .filter(ReportableDisruption::canonical)
                    .collect(Collectors.toList());
        }
        else
        {
            List<ReportableDisruption> disruptions = Lists.newArrayList();

            for (Disruption disruption : disruptionsOld)
            {
                if (reportableGenes.contains(disruption.gene()) && disruption.canonical() && disruption.isDisruptive())
                {
                    disruptions.add(ImmutableReportableDisruption.builder()
                            .svId(Integer.valueOf(disruption.svId()))
                            .chromosome(disruption.chromosome())
                            .orientation(disruption.orientation())
                            .strand(disruption.strand())
                            .chrBand(disruption.chrBand())
                            .gene(disruption.gene())
                            .canonical(disruption.canonical())
                            .type(disruption.type())
                            .ploidy(disruption.ploidy())
                            .exonUp(disruption.exonUp())
                            .exonDown(disruption.exonDown())
                            .build());
                }
            }

            return disruptions;
        }
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
