package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.ImmutableGeneDisruption;
import com.hartwig.hmftools.common.linx.ImmutableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.GeneDisruption;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.ImmutableGainLoss;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class CurationFunction {

    private CurationFunction() {
    }

    private static final String GENE_CDKN2A = "CDKN2A";
    private static final String GENE_CDKN2A_CANONICAL = "CDKN2A (p16)";
    private static final String GENE_CDKN2A_NON_CANONICAL = "CDKN2A (p14ARF)";

    @NotNull
    public static GenomicAnalysis curation(@NotNull GenomicAnalysis genomicAnalysis) {
        return ImmutableGenomicAnalysis.builder()
                .from(genomicAnalysis)
                .tumorSpecificEvidence(curateEvidence(genomicAnalysis.tumorSpecificEvidence()))
                .clinicalTrials(curateEvidence(genomicAnalysis.clinicalTrials()))
                .offLabelEvidence(curateEvidence(genomicAnalysis.offLabelEvidence()))
                .reportableVariants(curateReportableVariants(genomicAnalysis.reportableVariants()))
                .notifyGermlineStatusPerVariant(curateNotifyGermlineStatusPerVariant(genomicAnalysis.notifyGermlineStatusPerVariant()))
                .gainsAndLosses(curateGainsAndLosses(genomicAnalysis.gainsAndLosses()))
                .geneDisruptions(curateGeneDisruptions(genomicAnalysis.geneDisruptions()))
                .homozygousDisruptions(curateHomozygousDisruptions(genomicAnalysis.homozygousDisruptions()))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<ProtectEvidence> curateEvidence(@NotNull List<ProtectEvidence> tumorSpecificEvidences) {
        List<ProtectEvidence> curateTumorSpecificEvidences = Lists.newArrayList();
        for (ProtectEvidence evidence : tumorSpecificEvidences) {


            if (evidence.gene() != null && evidence.isCanonical() != null) {
                if (evidence.gene().equals(GENE_CDKN2A) && evidence.isCanonical()) {
                    curateTumorSpecificEvidences.add(ImmutableProtectEvidence.builder().from(evidence).gene(GENE_CDKN2A_CANONICAL).build());
                } else if (evidence.gene().equals(GENE_CDKN2A) && !evidence.isCanonical()) {
                    curateTumorSpecificEvidences.add(ImmutableProtectEvidence.builder()
                            .from(evidence)
                            .gene(GENE_CDKN2A_NON_CANONICAL)
                            .build());
                } else {
                    curateTumorSpecificEvidences.add(evidence);
                }
            } else {
                curateTumorSpecificEvidences.add(evidence);
            }

        }
        return curateTumorSpecificEvidences;
    }

    @NotNull
    @VisibleForTesting
    static List<ReportableVariant> curateReportableVariants(@NotNull List<ReportableVariant> reportableVariants) {
        List<ReportableVariant> curateReportableVariants = Lists.newArrayList();
        for (ReportableVariant variant : reportableVariants) {
            if (variant.gene().equals(GENE_CDKN2A) && variant.isCanonical()) {
                curateReportableVariants.add(ImmutableReportableVariant.builder().from(variant).gene(GENE_CDKN2A_CANONICAL).build());
            } else if (variant.gene().equals(GENE_CDKN2A) && !variant.isCanonical()) {
                curateReportableVariants.add(ImmutableReportableVariant.builder().from(variant).gene(GENE_CDKN2A_NON_CANONICAL).build());
            } else {
                curateReportableVariants.add(variant);
            }
        }

        return curateReportableVariants;
    }

    @NotNull
    @VisibleForTesting
    static Map<ReportableVariant, Boolean> curateNotifyGermlineStatusPerVariant(
            @NotNull Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariants) {
        Map<ReportableVariant, Boolean> curateNotifyGermlineStatusPerVariants = Maps.newHashMap();

        for (Map.Entry<ReportableVariant, Boolean> entry : notifyGermlineStatusPerVariants.entrySet()) {
            if (entry.getKey().gene().equals(GENE_CDKN2A) && entry.getKey().isCanonical()) {
                curateNotifyGermlineStatusPerVariants.put(ImmutableReportableVariant.builder()
                        .from(entry.getKey())
                        .gene(GENE_CDKN2A_CANONICAL)
                        .build(), entry.getValue());
            } else if (entry.getKey().gene().equals(GENE_CDKN2A) && !entry.getKey().isCanonical()) {
                curateNotifyGermlineStatusPerVariants.put(ImmutableReportableVariant.builder()
                        .from(entry.getKey())
                        .gene(GENE_CDKN2A_NON_CANONICAL)
                        .build(), entry.getValue());
            } else {
                curateNotifyGermlineStatusPerVariants.put(entry.getKey(), entry.getValue());
            }
        }

        return curateNotifyGermlineStatusPerVariants;
    }

    @NotNull
    @VisibleForTesting
    static List<GainLoss> curateGainsAndLosses(@NotNull List<GainLoss> gainsAndLosses) {
        List<GainLoss> curateGainsAndLosses = Lists.newArrayList();
        for (GainLoss gainLoss : gainsAndLosses) {
            if (gainLoss.gene().equals(GENE_CDKN2A) && gainLoss.isCanonical()) {
                curateGainsAndLosses.add(ImmutableGainLoss.builder().from(gainLoss).gene(GENE_CDKN2A_CANONICAL).build());
            } else if (gainLoss.gene().equals(GENE_CDKN2A) && !gainLoss.isCanonical()) {
                curateGainsAndLosses.add(ImmutableGainLoss.builder().from(gainLoss).gene(GENE_CDKN2A_NON_CANONICAL).build());
            } else {
                curateGainsAndLosses.add(gainLoss);
            }
        }

        return curateGainsAndLosses;
    }

    @NotNull
    @VisibleForTesting
    static List<GeneDisruption> curateGeneDisruptions(@NotNull List<GeneDisruption> geneDisruptions) {
        List<GeneDisruption> curateGeneDisruptions = Lists.newArrayList();
        for (GeneDisruption disruption : geneDisruptions) {
            if (disruption.gene().equals(GENE_CDKN2A) && disruption.isCanonical()) {
                curateGeneDisruptions.add(ImmutableGeneDisruption.builder().from(disruption).gene(GENE_CDKN2A_CANONICAL).build());
            } else if (disruption.gene().equals(GENE_CDKN2A) && !disruption.isCanonical()) {
                curateGeneDisruptions.add(ImmutableGeneDisruption.builder()
                        .from(disruption)
                        .gene(GENE_CDKN2A_NON_CANONICAL)
                        .build());
            } else {
                curateGeneDisruptions.add(disruption);
            }
        }
        return curateGeneDisruptions;
    }

    @NotNull
    @VisibleForTesting
    static List<HomozygousDisruption> curateHomozygousDisruptions(
            @NotNull List<HomozygousDisruption> homozygousDisruptions) {
        List<HomozygousDisruption> curateHomozygousDisruptions = Lists.newArrayList();
        for (HomozygousDisruption homozygousDisruption : homozygousDisruptions) {
            if (homozygousDisruption.gene().equals(GENE_CDKN2A) && homozygousDisruption.isCanonical()) {
                curateHomozygousDisruptions.add(ImmutableHomozygousDisruption.builder()
                        .from(homozygousDisruption)
                        .gene(GENE_CDKN2A_CANONICAL)
                        .build());
            } else if (homozygousDisruption.gene().equals(GENE_CDKN2A) && !homozygousDisruption.isCanonical()) {
                curateHomozygousDisruptions.add(ImmutableHomozygousDisruption.builder()
                        .from(homozygousDisruption)
                        .gene(GENE_CDKN2A_NON_CANONICAL)
                        .build());
            } else {
                curateHomozygousDisruptions.add(homozygousDisruption);
            }
        }
        return curateHomozygousDisruptions;
    }
}