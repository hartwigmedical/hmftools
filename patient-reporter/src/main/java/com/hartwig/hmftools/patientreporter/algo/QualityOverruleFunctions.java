package com.hartwig.hmftools.patientreporter.algo;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.patientreporter.QsFormNumber;

import org.jetbrains.annotations.NotNull;

public final class QualityOverruleFunctions {

    private QualityOverruleFunctions() {
    }

    @NotNull
    public static GenomicAnalysis overrule(@NotNull GenomicAnalysis genomicAnalysis, @NotNull String qcForm) {
        Map<ReportableVariant, Boolean> overruledVariantMaps =
                overruleMap(genomicAnalysis.notifyGermlineStatusPerVariant(), genomicAnalysis.hasReliablePurity());

        List<ReportableVariantWithNotify> overruledVariantsWithNotify =
                overruleVariants(genomicAnalysis.reportableVariants(), overruledVariantMaps, genomicAnalysis.hasReliablePurity());

        List<CnPerChromosomeArmData> cnPerChromosomeData = Lists.newArrayList();
        for (CnPerChromosomeArmData cnPerChromosome : genomicAnalysis.cnPerChromosome()) {
            if (genomicAnalysis.hasReliablePurity()) {
                cnPerChromosomeData.add(cnPerChromosome);
            }
        }

        List<ReportableVariant> overruledVariants = Lists.newArrayList();
        Map<ReportableVariant, Boolean> newNotifyPerVariant = Maps.newHashMap();
        for (ReportableVariantWithNotify overruled : overruledVariantsWithNotify) {
            overruledVariants.add(overruled.variant());
            newNotifyPerVariant.put(overruled.variant(), overruled.notifyVariant());
        }

        Map<ReportableVariant, Boolean> sortedNotifyPerVariant = newNotifyPerVariant.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByValue())
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        Map.Entry::getValue,
                        (oldValue, newValue) -> oldValue, LinkedHashMap::new));

        return ImmutableGenomicAnalysis.builder()
                .from(genomicAnalysis)
                .reportableVariants(overruledVariants)
                .notifyGermlineStatusPerVariant(sortedNotifyPerVariant)
                .cnPerChromosome(cnPerChromosomeData)
                .peachGenotypes(qcForm.equals(QsFormNumber.FOR_080.display()) ? genomicAnalysis.peachGenotypes() : Lists.newArrayList())
                .build();
    }

    @NotNull
    private static Map<ReportableVariant, Boolean> overruleMap(@NotNull Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant,
            boolean hasReliablePurity) {
        Map<ReportableVariant, Boolean> filteredMap = Maps.newHashMap();
        for (Map.Entry<ReportableVariant, Boolean> entry : notifyGermlineStatusPerVariant.entrySet()) {
            filteredMap.put(ImmutableReportableVariant.builder()
                    .from(QualityOverruleFunctions.overruleVariant(entry.getKey(), hasReliablePurity))
                    .build(), entry.getValue());

        }
        return filteredMap;
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
                .minorAlleleCopyNumber(hasReliablePurity ? variant.minorAlleleCopyNumber() : Double.NaN)
                .biallelic(hasReliablePurity && roundedCopyNumber >= 1 ? variant.biallelic() : null)
                .build();
    }
}
