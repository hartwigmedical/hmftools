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
        List<ReportableVariantWithNotify> overruledVariantsWithNotify = overruleVariants(genomicAnalysis.reportableVariants(),
                genomicAnalysis.notifyGermlineStatusPerVariant(),
                genomicAnalysis.hasReliablePurity());

        List<ReportableVariant> overruledVariants = Lists.newArrayList();
        Map<ReportableVariant, Boolean> newNotifyPerVariant = Maps.newHashMap();
        for (ReportableVariantWithNotify overruled : overruledVariantsWithNotify) {
            overruledVariants.add(overruled.variant());
            newNotifyPerVariant.put(overruled.variant(), overruled.notifyVariant());
        }

        return ImmutableGenomicAnalysis.builder()
                .from(genomicAnalysis)
                .reportableVariants(overruledVariants)
                .notifyGermlineStatusPerVariant(newNotifyPerVariant)
                .build();
    }

    @NotNull
    private static List<ReportableVariantWithNotify> overruleVariants(@NotNull List<ReportableVariant> variants,
            @NotNull Map<ReportableVariant, Boolean> notifyVariant, boolean hasReliablePurity) {
        List<ReportableVariantWithNotify> overruledVariants = Lists.newArrayList();

        for (ReportableVariant variant : variants) {
            ReportableVariant newVariant = overruleVariant(variant, hasReliablePurity);
            overruledVariants.add(ImmutableReportableVariantWithNotify.builder()
                    .variant(newVariant)
                    .notifyVariant(notifyVariant.get(newVariant))
                    .build());
        }

        return overruledVariants;
    }

    @NotNull
    private static ReportableVariant overruleVariant(@NotNull ReportableVariant variant, boolean hasReliablePurity) {
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
