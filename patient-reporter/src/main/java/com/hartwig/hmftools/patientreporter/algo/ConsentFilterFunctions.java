package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;

import org.immutables.value.internal.$guava$.annotations.$VisibleForTesting;
import org.jetbrains.annotations.NotNull;

public final class ConsentFilterFunctions {

    private ConsentFilterFunctions() {
    }

    @NotNull
    public static GenomicAnalysis filter(@NotNull GenomicAnalysis genomicAnalysis,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel, boolean reportViralPresence, boolean reportPeach) {
        List<ReportableVariantWithNotify> filteredVariantsWithNotify = filterVariants(genomicAnalysis.reportableVariants(),
                genomicAnalysis.notifyGermlineStatusPerVariant(),
                germlineReportingLevel);

        List<ReportableVariant> filteredVariants = Lists.newArrayList();
        Map<ReportableVariant, Boolean> notifyPerVariant = Maps.newHashMap();
        for (ReportableVariantWithNotify filtered : filteredVariantsWithNotify) {
            filteredVariants.add(filtered.variant());
            notifyPerVariant.put(filtered.variant(), filtered.notifyVariant());
        }

        List<AnnotatedVirus> filteredViruses =
                reportViralPresence ? genomicAnalysis.reportableViruses() : Lists.newArrayList();

        List<PeachGenotype> filteredPeachGenotypes = reportPeach ? genomicAnalysis.peachGenotypes() : Lists.newArrayList();

        List<ProtectEvidence> filteredTumorSpecificEvidence =
                filterEvidenceForGermlineConsent(genomicAnalysis.tumorSpecificEvidence(), germlineReportingLevel);

        List<ProtectEvidence> filteredClinicalTrials =
                filterEvidenceForGermlineConsent(genomicAnalysis.clinicalTrials(), germlineReportingLevel);

        List<ProtectEvidence> filteredOffLabelEvidence =
                filterEvidenceForGermlineConsent(genomicAnalysis.offLabelEvidence(), germlineReportingLevel);

        return ImmutableGenomicAnalysis.builder()
                .from(genomicAnalysis)
                .reportableVariants(filteredVariants)
                .notifyGermlineStatusPerVariant(notifyPerVariant)
                .reportableViruses(filteredViruses)
                .peachGenotypes(filteredPeachGenotypes)
                .tumorSpecificEvidence(filteredTumorSpecificEvidence)
                .clinicalTrials(filteredClinicalTrials)
                .offLabelEvidence(filteredOffLabelEvidence)
                .build();
    }

    @NotNull
    @$VisibleForTesting
    static List<ReportableVariantWithNotify> filterVariants(@NotNull List<ReportableVariant> variants,
            @NotNull Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        List<ReportableVariantWithNotify> filteredVariants = Lists.newArrayList();
        for (ReportableVariant variant : variants) {
            if (!(variant.source() == ReportableVariantSource.GERMLINE
                    && germlineReportingLevel == LimsGermlineReportingLevel.NO_REPORTING)) {
                if (variant.source() == ReportableVariantSource.GERMLINE && !notifyGermlineStatusPerVariant.get(variant)) {
                    filteredVariants.add(ImmutableReportableVariantWithNotify.builder()
                            .variant(ImmutableReportableVariant.builder().from(variant).source(ReportableVariantSource.SOMATIC).build())
                            .notifyVariant(false)
                            .build());
                } else {
                    filteredVariants.add(ImmutableReportableVariantWithNotify.builder()
                            .variant(variant)
                            .notifyVariant(notifyGermlineStatusPerVariant.get(variant))
                            .build());
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
            if (evidence.germline() && germlineReportingLevel != LimsGermlineReportingLevel.NO_REPORTING) {
                // We always overwrite to somatic in evidence since we are not sure we notify about the actual variant.
                filtered.add(ImmutableProtectEvidence.builder().from(evidence).germline(false).build());
            } else if (!evidence.germline()) {
                filtered.add(evidence);
            }
        }

        return filtered;
    }
}
