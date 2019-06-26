package com.hartwig.hmftools.patientreporter.structural;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SvAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(SvAnalyzer.class);

    @NotNull
    private final List<ReportableGeneFusion> fusions;
    @NotNull
    private final List<ReportableDisruption> disruptions;

    @NotNull
    public static SvAnalyzer fromFiles(@NotNull String linxFusionTsv, @NotNull String linxDisruptionTsv) throws IOException {
        LOGGER.info("Loading fusions from {}", linxFusionTsv);
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.read(linxFusionTsv);
        LOGGER.info(" Loaded {} fusions", fusions.size());

        LOGGER.info("Loading disruptions from {}", linxDisruptionTsv);
        List<ReportableDisruption> disruptions = ReportableDisruptionFile.read(linxDisruptionTsv);
        LOGGER.info(" Loaded {} disruptions", disruptions.size());

        return new SvAnalyzer(fusions, disruptions);
    }

    @VisibleForTesting
    SvAnalyzer(@NotNull final List<ReportableGeneFusion> fusions, @NotNull final List<ReportableDisruption> disruptions) {
        this.fusions = fusions;
        this.disruptions = disruptions;
    }

    @NotNull
    public SvAnalysis run(@NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation) {
        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(disruptions, geneCopyNumbers);

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;

        Map<SimpleGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(toSimpleGeneFusions(fusions), primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerFusion);

        return ImmutableSvAnalysis.builder()
                .reportableFusions(fusions)
                .reportableDisruptions(reportableGeneDisruptions)
                .evidenceItems(filteredEvidence)
                .build();
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
