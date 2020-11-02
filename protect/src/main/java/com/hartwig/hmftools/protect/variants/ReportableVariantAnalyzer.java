package com.hartwig.hmftools.protect.variants;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.protect.actionability.ReportableEvidenceItemFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableVariantAnalyzer {

    private ReportableVariantAnalyzer() {
    }

    @NotNull
    public static ReportableVariantAnalysis mergeSomaticAndGermlineVariants(@NotNull List<ReportableVariant> reportableSomaticVariants,
            @NotNull List<ReportableVariant> reportableGermlineVariants, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation) {
        List<ReportableVariant> allReportableVariants =
                ReportableVariantFactory.mergeVariantLists(reportableGermlineVariants, reportableSomaticVariants);

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        // Extract somatic evidence for high drivers variants only (See DEV-824)
        Map<ReportableVariant, List<EvidenceItem>> evidencePerVariant =
                filterHighDriverLikelihood(actionabilityAnalyzer.evidenceForAllVariants(allReportableVariants, primaryTumorLocation));

        return ImmutableReportableVariantAnalysis.of(allReportableVariants,
                ReportableEvidenceItemFactory.toReportableFlatList(evidencePerVariant));
    }

    @NotNull
    private static Map<ReportableVariant, List<EvidenceItem>> filterHighDriverLikelihood(
            @NotNull Map<? extends Variant, List<EvidenceItem>> evidenceForAllVariants) {
        Map<ReportableVariant, List<EvidenceItem>> evidencePerHighDriverVariant = Maps.newHashMap();
        for (Map.Entry<? extends Variant, List<EvidenceItem>> entry : evidenceForAllVariants.entrySet()) {
            ReportableVariant variant = (ReportableVariant) entry.getKey();
            if (DriverInterpretation.interpret(variant.driverLikelihood()) == DriverInterpretation.HIGH) {
                evidencePerHighDriverVariant.put(variant, entry.getValue());
            }
        }
        return evidencePerHighDriverVariant;
    }
}

