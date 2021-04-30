package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.protect.viralbreakend.VirusBreakend;
import com.hartwig.hmftools.protect.purple.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.immutables.value.internal.$guava$.annotations.$VisibleForTesting;
import org.jetbrains.annotations.NotNull;

public final class ConsentFilterFunctions {

    public ConsentFilterFunctions() {
    }

    // TODO Split up the filtering functions from the overruling functions.
    // TODO Extend overruling functions to make json and pdf more consistent.
    @NotNull
    public GenomicAnalysis filterAndOverruleForConsent(@NotNull GenomicAnalysis genomicAnalysis,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel, boolean reportViralBreakends) {
        List<ReportableVariant> filteredVariants = filterAndOverruleVariants(genomicAnalysis.reportableVariants(),
                germlineReportingLevel,
                genomicAnalysis.hasReliablePurity());

        List<VirusBreakend> filteredVirusBreakends = reportViralBreakends ? genomicAnalysis.virusBreakends() : Lists.newArrayList();

        List<ProtectEvidence> filteredTumorSpecificEvidence =
                filterEvidenceForGermlineConsent(genomicAnalysis.tumorSpecificEvidence(), germlineReportingLevel);

        List<ProtectEvidence> filteredClinicalTrials =
                filterEvidenceForGermlineConsent(genomicAnalysis.clinicalTrials(), germlineReportingLevel);

        List<ProtectEvidence> filteredOffLabelEvidence =
                filterEvidenceForGermlineConsent(genomicAnalysis.offLabelEvidence(), germlineReportingLevel);

        return ImmutableGenomicAnalysis.builder()
                .from(genomicAnalysis)
                .reportableVariants(filteredVariants)
                .virusBreakends(filteredVirusBreakends)
                .tumorSpecificEvidence(filteredTumorSpecificEvidence)
                .clinicalTrials(filteredClinicalTrials)
                .offLabelEvidence(filteredOffLabelEvidence)
                .build();
    }

    @NotNull
    @$VisibleForTesting
    static List<ReportableVariant> filterAndOverruleVariants(@NotNull List<ReportableVariant> variants,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel, boolean hasReliablePurity) {

        List<ReportableVariant> filteredVariants = Lists.newArrayList();
        for (ReportableVariant variant : variants) {
            if (germlineReportingLevel != LimsGermlineReportingLevel.NO_REPORTING || variant.source() == ReportableVariantSource.SOMATIC) {
                filteredVariants.add(variant);
            }
        }

        List<ReportableVariant> overruledVariants = Lists.newArrayList();
        for (ReportableVariant variant : filteredVariants) {
            if (germlineReportingLevel == LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION) {
                overruledVariants.add(overruleVariant(variant, hasReliablePurity, ReportableVariantSource.SOMATIC));
            } else {
                overruledVariants.add(overruleVariant(variant, hasReliablePurity, variant.source()));
            }
        }
        return overruledVariants;
    }

    @NotNull
    @VisibleForTesting
    static List<ProtectEvidence> filterEvidenceForGermlineConsent(@NotNull List<ProtectEvidence> evidences,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.germline() && germlineReportingLevel == LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION) {
                filtered.add(ImmutableProtectEvidence.builder().from(evidence).germline(false).build());
            } else if (evidence.germline() && germlineReportingLevel == LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION) {
                filtered.add(evidence);
            } else if (!evidence.germline()) {
                filtered.add(evidence);
            }
        }

        return filtered;
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
                .biallelic(hasReliablePurity ? variant.biallelic() : null)
                .build();
    }
}
