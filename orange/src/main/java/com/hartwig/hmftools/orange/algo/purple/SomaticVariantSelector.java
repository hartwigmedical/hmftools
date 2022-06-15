package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.orange.algo.protect.EvidenceSelector;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class SomaticVariantSelector {

    static final Set<String> CUPPA_GENES = Sets.newHashSet("ALB", "SFTPB", "SLC34A2");

    private SomaticVariantSelector() {
    }

    @NotNull
    public static List<ReportableVariant> selectInterestingUnreportedVariants(@NotNull List<SomaticVariant> allVariants,
            @NotNull List<ReportableVariant> reportedSomaticVariants, @NotNull List<ProtectEvidence> evidences,
            @NotNull List<DriverGene> driverGenes) {
        List<ReportableVariant> filtered = Lists.newArrayList();
        for (SomaticVariant variant : allVariants) {
            if (!variant.reported()) {
                boolean isNearHotspot = variant.hotspot() == Hotspot.HOTSPOT || variant.hotspot() == Hotspot.NEAR_HOTSPOT;
                boolean hasEvidence = EvidenceSelector.hasEvidence(evidences, variant.gene(), ProtectEventGenerator.variantEvent(variant));
                boolean isExonicAndHasPhasedReportedVariant =
                        !variant.gene().isEmpty() && hasReportedVariantWithPhase(reportedSomaticVariants, variant.topLocalPhaseSet());
                boolean isCuppaRelevantVariant = isRelevantForCuppa(variant);
                boolean isSynonymousButReportable = isSynonymousWithReportableWorstImpact(variant, driverGenes);
                boolean isUnreportedSpliceVariant = isUnreportedSpliceVariant(variant, driverGenes);

                if (isNearHotspot || hasEvidence || isExonicAndHasPhasedReportedVariant || isCuppaRelevantVariant
                        || isSynonymousButReportable || isUnreportedSpliceVariant) {
                    filtered.add(toReportable(variant));
                }
            }
        }
        return filtered;
    }

    private static boolean hasReportedVariantWithPhase(@NotNull List<ReportableVariant> reportedVariants,
            @Nullable Integer targetPhaseSet) {
        if (targetPhaseSet == null) {
            return false;
        }

        for (ReportableVariant variant : reportedVariants) {
            if (variant.localPhaseSet() != null && variant.localPhaseSet().equals(targetPhaseSet)) {
                return true;
            }
        }

        return false;
    }

    private static boolean isRelevantForCuppa(@NotNull SomaticVariant variant) {
        return variant.type() == VariantType.INDEL && CUPPA_GENES.contains(variant.gene()) && variant.repeatCount() <= 6;
    }

    private static boolean isSynonymousWithReportableWorstImpact(@NotNull SomaticVariant variant, @NotNull List<DriverGene> driverGenes) {
        if (variant.canonicalCodingEffect() != CodingEffect.SYNONYMOUS) {
            return false;
        }

        DriverGene driverGene = findDriverGene(driverGenes, variant.gene());
        if (driverGene == null) {
            return false;
        }

        CodingEffect worstEffect = variant.worstCodingEffect();
        boolean nonsenseOrFrameshift = worstEffect == CodingEffect.NONSENSE_OR_FRAMESHIFT && driverGene.reportNonsenseAndFrameshift();
        boolean splice = worstEffect == CodingEffect.SPLICE && driverGene.reportSplice();
        boolean missense = worstEffect == CodingEffect.MISSENSE && driverGene.reportMissenseAndInframe();
        return nonsenseOrFrameshift || splice || missense;
    }

    private static boolean isUnreportedSpliceVariant(@NotNull SomaticVariant variant, @NotNull List<DriverGene> driverGenes) {
        if (variant.spliceRegion()) {
            DriverGene driverGene = findDriverGene(driverGenes, variant.gene());
            if (driverGene != null) {
                return driverGene.reportSplice();
            }
        }
        return false;
    }

    @Nullable
    private static DriverGene findDriverGene(@NotNull List<DriverGene> driverGenes, @NotNull String geneToFind) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(geneToFind)) {
                return driverGene;
            }
        }
        return null;
    }

    @NotNull
    private static ReportableVariant toReportable(@NotNull SomaticVariant variant) {
        return ReportableVariantFactory.fromVariant(variant, ReportableVariantSource.SOMATIC)
                .driverLikelihood(Double.NaN)
                .transcript(variant.canonicalTranscript())
                .isCanonical(true)
                .build();
    }
}
