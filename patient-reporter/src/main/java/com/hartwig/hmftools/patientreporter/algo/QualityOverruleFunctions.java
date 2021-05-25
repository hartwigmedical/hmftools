package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.protect.purple.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariant;

import org.jetbrains.annotations.NotNull;

public final class QualityOverruleFunctions {

    private QualityOverruleFunctions() {
    }

    @NotNull
    public static GenomicAnalysis overrule(@NotNull GenomicAnalysis genomicAnalysis) {
        List<ReportableVariant> overruledVariants =
                overruleVariants(genomicAnalysis.reportableVariants(), genomicAnalysis.hasReliablePurity());

        return ImmutableGenomicAnalysis.builder()
                .from(genomicAnalysis)
                .reportableVariants(overruledVariants)
                .notifyGermlineStatusPerVariant(filterVariantsNotifyMap(genomicAnalysis.notifyGermlineStatusPerVariant(),
                        genomicAnalysis.hasReliablePurity())).build();
    }

    private static Map<ReportableVariant, Boolean> filterVariantsNotifyMap(
            @NotNull Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant, boolean hasReliablePurity) {

        Map<ReportableVariant, Boolean> filteredMap = Maps.newHashMap();
        for (Map.Entry<ReportableVariant, Boolean> entry : notifyGermlineStatusPerVariant.entrySet()) {
            filteredMap.put(ImmutableReportableVariant.builder()
                    .from(QualityOverruleFunctions.overruleVariant(entry.getKey(), hasReliablePurity))
                    .source(entry.getKey().source())
                    .build(), entry.getValue());
        }
        return filteredMap;
    }

    @NotNull
    private static List<ReportableVariant> overruleVariants(@NotNull List<ReportableVariant> variants, boolean hasReliablePurity) {
        List<ReportableVariant> overruledVariants = Lists.newArrayList();

        for (ReportableVariant variant : variants) {
            overruledVariants.add(overruleVariant(variant, hasReliablePurity));
        }

        return overruledVariants;
    }

    @NotNull
    public static ReportableVariant overruleVariant(@NotNull ReportableVariant variant, boolean hasReliablePurity) {
        double flooredCopyNumber = Math.max(0, variant.totalCopyNumber());
        long roundedCopyNumber = Math.round(flooredCopyNumber);

        return ImmutableReportableVariant.builder()
                .from(variant)
                .totalCopyNumber(hasReliablePurity && roundedCopyNumber >= 1 ? flooredCopyNumber : Double.NaN)
                .alleleCopyNumber(hasReliablePurity && roundedCopyNumber >= 1 ? variant.alleleCopyNumber() : Double.NaN)
                .biallelic(hasReliablePurity && roundedCopyNumber >= 1 ? variant.biallelic() : null)
                .build();
    }
}
