package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.reportablegenomicalterations.AllReportableVariants;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.common.drivergene.DriverGeneView;
import com.hartwig.hmftools.common.bachelor.GermlineReportingModel;
import com.hartwig.hmftools.common.reportablegenomicalterations.ReportableGermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableVariantAnalyzer {

    private ReportableVariantAnalyzer() {
    }

    @NotNull
    public static ReportVariantAnalysis mergeSomaticAndGermlineVariants(@NotNull List<SomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull DriverGeneView driverGeneView,
            @NotNull List<ReportableGermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingChoice germlineReportingChoice, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation) {
        List<ReportableVariant> allReportableVariants = AllReportableVariants.mergeSomaticAndGermlineVariants(somaticVariantsReport,
                driverCatalog,
                driverGeneView,
                germlineVariantsToReport,
                germlineReportingModel, germlineReportingChoice);

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        // Extract somatic evidence for high drivers variants only (See DEV-824)
        Map<ReportableVariant, List<EvidenceItem>> evidencePerVariant =
                filterHighDriverLikelihood(actionabilityAnalyzer.evidenceForAllVariants(allReportableVariants, primaryTumorLocation));

        return ImmutableReportVariantAnalysis.of(allReportableVariants,
                ReportableEvidenceItemFactory.toReportableFlatList(evidencePerVariant));
    }

    @NotNull
    private static Map<ReportableVariant, List<EvidenceItem>> filterHighDriverLikelihood(
            final Map<ReportableVariant, List<EvidenceItem>> evidenceForAllVariants) {
        Map<ReportableVariant, List<EvidenceItem>> evidencePerHighDriverVariant = Maps.newHashMap();
        for (Map.Entry<ReportableVariant, List<EvidenceItem>> entry : evidenceForAllVariants.entrySet()) {
            if (DriverInterpretation.interpret(entry.getKey().driverLikelihood()) == DriverInterpretation.HIGH) {
                evidencePerHighDriverVariant.put(entry.getKey(), entry.getValue());
            }
        }
        return evidencePerHighDriverVariant;
    }
}

