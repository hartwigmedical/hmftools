package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.protect.purple.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.jetbrains.annotations.NotNull;

public final class ConsentFilterFunctions {

    private ConsentFilterFunctions() {
    }

    @NotNull
    public static AnalysedPatientReport filterForConsent(@NotNull AnalysedPatientReport report) {
        List<ReportableVariant> filteredVariantsOverruleVariantSource = Lists.newArrayList();
        for (ReportableVariant variant : report.genomicAnalysis().reportableVariants()) {
            if (report.sampleReport().germlineReportingLevel() == LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION) {
                filteredVariantsOverruleVariantSource.add(overruleVariant(variant,
                        report.genomicAnalysis().hasReliablePurity(),
                        ReportableVariantSource.SOMATIC));
            } else {
                filteredVariantsOverruleVariantSource.add(overruleVariant(variant,
                        report.genomicAnalysis().hasReliablePurity(),
                        variant.source()));
            }
        }

        return ImmutableAnalysedPatientReport.builder()
                .from(report)
                .genomicAnalysis(replaceReportableVariants(report, filteredVariantsOverruleVariantSource))
                .build();
    }

    @NotNull
    private static ImmutableGenomicAnalysis replaceReportableVariants(@NotNull AnalysedPatientReport report,
            @NotNull List<ReportableVariant> filteredVariantsOverruleVariantSource) {
        return ImmutableGenomicAnalysis.builder()
                .from(report.genomicAnalysis())
                .reportableVariants(filteredVariantsOverruleVariantSource)
                .build();
    }

    @NotNull
    private static ReportableVariant overruleVariant(@NotNull ReportableVariant variant, boolean hasReliablePurity,
            @NotNull ReportableVariantSource source) {
        double flooredCopyNumber = Math.max(0, variant.totalCopyNumber());

        return ImmutableReportableVariant.builder()
                .from(variant)
                .source(source)
                .totalCopyNumber(hasReliablePurity ? flooredCopyNumber : Double.NaN)
                .alleleCopyNumber(hasReliablePurity ? variant.alleleCopyNumber() : Double.NaN)
                .build();
    }
}
