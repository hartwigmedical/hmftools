package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakend;
import com.hartwig.hmftools.protect.purple.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.immutables.value.internal.$guava$.annotations.$VisibleForTesting;
import org.jetbrains.annotations.NotNull;

public final class ConsentFilterFunctions {

    private ConsentFilterFunctions() {
    }

    @NotNull
    public static GenomicAnalysis filter(@NotNull GenomicAnalysis genomicAnalysis,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel, boolean reportViralBreakends, boolean reportPeach) {
        List<ReportableVariant> filteredVariants = filterVariants(genomicAnalysis.reportableVariants(), germlineReportingLevel);

        List<ReportableVirusBreakend> filteredVirusBreakends =
                reportViralBreakends ? genomicAnalysis.virusBreakends() : Lists.newArrayList();

        List<PeachGenotype> filteredPeachGenotype = reportPeach ? genomicAnalysis.peachGenotypes() : Lists.newArrayList();

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
                .peachGenotypes(filteredPeachGenotype)
                .tumorSpecificEvidence(filteredTumorSpecificEvidence)
                .clinicalTrials(filteredClinicalTrials)
                .offLabelEvidence(filteredOffLabelEvidence)
                .build();
    }

    @NotNull
    @$VisibleForTesting
    static List<ReportableVariant> filterVariants(@NotNull List<ReportableVariant> variants,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        List<ReportableVariant> filteredVariants = Lists.newArrayList();
        for (ReportableVariant variant : variants) {
            if (germlineReportingLevel != LimsGermlineReportingLevel.NO_REPORTING || variant.source() == ReportableVariantSource.SOMATIC) {
                if (germlineReportingLevel == LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION
                        && variant.source() == ReportableVariantSource.GERMLINE) {
                    filteredVariants.add(ImmutableReportableVariant.builder()
                            .from(variant)
                            .source(ReportableVariantSource.SOMATIC)
                            .build());
                } else {
                    filteredVariants.add(variant);
                }
            }
        }
        return filteredVariants;
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
}
