package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class SomaticVariantSelector {

    static final Set<String> CUPPA_GENES = Sets.newHashSet("ALB", "SFTPB", "SLC34A2");

    private SomaticVariantSelector() {
    }

    @NotNull
    public static List<PurpleVariant> selectInterestingUnreportedVariants(@NotNull List<PurpleVariant> allSomaticVariants,
            @NotNull List<PurpleVariant> reportedSomaticVariants, @NotNull List<DriverGene> driverGenes) {
        List<PurpleVariant> filtered = Lists.newArrayList();
        for (PurpleVariant variant : allSomaticVariants) {
            if (!variant.reported()) {
                boolean isNearHotspot = variant.hotspot() == Hotspot.HOTSPOT || variant.hotspot() == Hotspot.NEAR_HOTSPOT;
                boolean isExonicAndHasPhasedReportedVariant =
                        !variant.gene().isEmpty() && hasReportedVariantWithPhase(reportedSomaticVariants, variant.localPhaseSets());
                boolean isCuppaRelevantVariant = isRelevantForCuppa(variant);
                boolean isSynonymousButReportable = isSynonymousWithReportableWorstImpact(variant, driverGenes);
                boolean isUnreportedSpliceVariant = isUnreportedSpliceVariant(variant, driverGenes);

                if (isNearHotspot || isExonicAndHasPhasedReportedVariant || isCuppaRelevantVariant || isSynonymousButReportable
                        || isUnreportedSpliceVariant) {
                    filtered.add(variant);
                }
            }
        }
        return filtered;
    }

    private static boolean hasReportedVariantWithPhase(@NotNull List<PurpleVariant> reportedVariants,
            @Nullable List<Integer> targetPhaseSets) {
        if (targetPhaseSets == null) {
            return false;
        }

        for (PurpleVariant variant : reportedVariants) {
            List<Integer> localPhaseSets = variant.localPhaseSets();
            if (localPhaseSets != null && hasMatchingPhase(localPhaseSets, targetPhaseSets)) {
                return true;
            }
        }

        return false;
    }

    private static boolean hasMatchingPhase(@NotNull List<Integer> localPhaseSets, @NotNull List<Integer> targetPhaseSets) {
        for (Integer localPhaseSet : localPhaseSets) {
            if (targetPhaseSets.contains(localPhaseSet)) {
                return true;
            }
        }
        return false;
    }

    private static boolean isRelevantForCuppa(@NotNull PurpleVariant variant) {
        return variant.type() == VariantType.INDEL && CUPPA_GENES.contains(variant.gene()) && variant.repeatCount() <= 6;
    }

    private static boolean isSynonymousWithReportableWorstImpact(@NotNull PurpleVariant variant, @NotNull List<DriverGene> driverGenes) {
        if (variant.canonicalImpact().codingEffect() != CodingEffect.SYNONYMOUS) {
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

    private static boolean isUnreportedSpliceVariant(@NotNull PurpleVariant variant, @NotNull List<DriverGene> driverGenes) {
        if (variant.canonicalImpact().spliceRegion()) {
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
}
