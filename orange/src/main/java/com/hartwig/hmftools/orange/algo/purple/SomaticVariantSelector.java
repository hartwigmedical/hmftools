package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.cuppa.CuppaCommon;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class SomaticVariantSelector
{
    static final Set<String> CUPPA_GENES = Sets.newHashSet(
            CuppaCommon.INDEL_ALB, CuppaCommon.INDEL_SFTPB, CuppaCommon.INDEL_SLC34A2);

    public static List<PurpleVariant> selectInterestingUnreportedVariants(
            final List<PurpleVariant> allSomaticVariants, final List<PurpleVariant> reportedSomaticVariants, final List<DriverGene> driverGenes)
    {
        List<PurpleVariant> filtered = Lists.newArrayList();
        for(PurpleVariant variant : allSomaticVariants)
        {
            if(!variant.reported())
            {
                boolean isAtLeastNearHotspot = variant.hotspot() == HotspotType.HOTSPOT || variant.hotspot() == HotspotType.NEAR_HOTSPOT;
                boolean isExonicAndHasPhasedReportedVariant =
                        !variant.gene().isEmpty() && isExonic(variant) && hasReportedVariantWithPhase(reportedSomaticVariants, variant.localPhaseSets());
                boolean isCuppaRelevantVariant = isRelevantForCuppa(variant);
                boolean isSynonymousButReportable = isSynonymousWithReportableWorstImpact(variant, driverGenes);
                boolean isUnreportedSpliceVariant = isUnreportedSpliceVariant(variant, driverGenes);

                if(isAtLeastNearHotspot || isExonicAndHasPhasedReportedVariant || isCuppaRelevantVariant || isSynonymousButReportable
                || isUnreportedSpliceVariant)
                {
                    filtered.add(variant);
                }
            }
        }
        return filtered;
    }

    private static boolean hasReportedVariantWithPhase(final List<PurpleVariant> reportedVariants, @Nullable List<Integer> targetPhaseSets)
    {
        if(targetPhaseSets == null)
        {
            return false;
        }

        for(PurpleVariant variant : reportedVariants)
        {
            List<Integer> localPhaseSets = variant.localPhaseSets();
            if(localPhaseSets != null && hasMatchingPhase(localPhaseSets, targetPhaseSets))
            {
                return true;
            }
        }

        return false;
    }

    private static boolean hasMatchingPhase(final List<Integer> localPhaseSets, final List<Integer> targetPhaseSets)
    {
        for(Integer localPhaseSet : localPhaseSets)
        {
            if(targetPhaseSets.contains(localPhaseSet))
            {
                return true;
            }
        }
        return false;
    }

    private static boolean isExonic(final PurpleVariant variant) {
        return variant.canonicalImpact().affectedExon() != null;
    }

    private static boolean isRelevantForCuppa(final PurpleVariant variant)
    {
        return variant.type() == PurpleVariantType.INDEL
            && CUPPA_GENES.contains(variant.gene())
            && variant.repeatCount() <= CuppaCommon.INDEL_MAX_REPEAT_COUNT;
    }

    private static boolean isSynonymousWithReportableWorstImpact(final PurpleVariant variant, final List<DriverGene> driverGenes)
    {
        if(variant.canonicalImpact().codingEffect() != PurpleCodingEffect.SYNONYMOUS)
        {
            return false;
        }

        DriverGene driverGene = findDriverGene(driverGenes, variant.gene());
        if(driverGene == null)
        {
            return false;
        }

        PurpleCodingEffect worstEffect = variant.worstCodingEffect();
        boolean nonsenseOrFrameshift = worstEffect == PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT && driverGene.reportNonsenseAndFrameshift();
        boolean splice = worstEffect == PurpleCodingEffect.SPLICE && driverGene.reportSplice();
        boolean missense = worstEffect == PurpleCodingEffect.MISSENSE && driverGene.reportMissenseAndInframe();
        return nonsenseOrFrameshift || splice || missense;
    }

    private static boolean isUnreportedSpliceVariant(final PurpleVariant variant, final List<DriverGene> driverGenes)
    {
        if(variant.canonicalImpact().inSpliceRegion())
        {
            DriverGene driverGene = findDriverGene(driverGenes, variant.gene());
            if(driverGene != null)
            {
                return driverGene.reportSplice();
            }
        }
        return false;
    }

    @Nullable
    private static DriverGene findDriverGene(final List<DriverGene> driverGenes, @NotNull String geneToFind)
    {
        for(DriverGene driverGene : driverGenes)
        {
            if(driverGene.gene().equals(geneToFind))
            {
                return driverGene;
            }
        }
        return null;
    }
}
