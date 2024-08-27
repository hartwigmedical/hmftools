package com.hartwig.hmftools.orange.algo.util;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportLimiter
{
    @NotNull
    public static OrangeRecord limitAllListsToMaxOne(@NotNull OrangeRecord report)
    {
        return ImmutableOrangeRecord.builder()
                .from(report)
                .germlineMVLHPerGene(limitGermlineMVLHToOne(report.germlineMVLHPerGene()))
                .purple(limitPurpleDataToOne(report.purple()))
                .linx(limitLinxDataToOne(report.linx()))
                .isofox(limitIsofoxDataToOne(report.isofox()))
                .build();
    }

    @Nullable
    private static Map<String, Double> limitGermlineMVLHToOne(@Nullable Map<String, Double> germlineMVLHPerGene)
    {
        if(germlineMVLHPerGene == null)
        {
            return null;
        }

        Map<String, Double> filtered = Maps.newHashMap();
        if(!germlineMVLHPerGene.isEmpty())
        {
            String firstKey = germlineMVLHPerGene.keySet().iterator().next();
            filtered.put(firstKey, germlineMVLHPerGene.get(firstKey));
        }
        return filtered;
    }

    @NotNull
    private static PurpleRecord limitPurpleDataToOne(@NotNull PurpleRecord purple)
    {
        return ImmutablePurpleRecord.builder()
                .from(purple)
                .somaticDrivers(max1(purple.somaticDrivers()))
                .germlineDrivers(max1(purple.germlineDrivers()))
                .allSomaticVariants(max1(purple.allSomaticVariants()))
                .reportableSomaticVariants(max1(purple.reportableSomaticVariants()))
                .additionalSuspectSomaticVariants(max1(purple.additionalSuspectSomaticVariants()))
                .allGermlineVariants(max1(purple.allGermlineVariants()))
                .reportableGermlineVariants(max1(purple.reportableGermlineVariants()))
                .additionalSuspectGermlineVariants(max1(purple.additionalSuspectGermlineVariants()))
                .allSomaticCopyNumbers(max1(purple.allSomaticCopyNumbers()))
                .allSomaticGeneCopyNumbers(max1(purple.allSomaticGeneCopyNumbers()))
                .suspectGeneCopyNumbersWithLOH(max1(purple.suspectGeneCopyNumbersWithLOH()))
                .allSomaticGainsLosses(max1(purple.allSomaticGainsLosses()))
                .reportableSomaticGainsLosses(max1(purple.reportableSomaticGainsLosses()))
                .nearReportableSomaticGains(max1(purple.nearReportableSomaticGains()))
                .additionalSuspectSomaticGainsLosses(max1(purple.additionalSuspectSomaticGainsLosses()))
                .allGermlineDeletions(max1(purple.allGermlineDeletions()))
                .allGermlineFullLosses(max1(purple.allGermlineFullLosses()))
                .reportableGermlineFullLosses(max1(purple.reportableGermlineFullLosses()))
                .allGermlineLossOfHeterozygosities(max1(purple.allGermlineLossOfHeterozygosities()))
                .reportableGermlineLossOfHeterozygosities(max1(purple.reportableGermlineLossOfHeterozygosities()))
                .build();
    }

    @NotNull
    private static LinxRecord limitLinxDataToOne(@NotNull LinxRecord linx)
    {
        List<LinxBreakend> filteredAllSomaticBreakends = max1(linx.allSomaticBreakends());
        List<LinxBreakend> filteredReportableSomaticBreakends = max1(linx.reportableSomaticBreakends());
        List<LinxBreakend> filteredAllGermlineBreakends = max1(linx.allGermlineBreakends());
        List<LinxBreakend> filteredReportableGermlineBreakends = max1(linx.reportableGermlineBreakends());

        return ImmutableLinxRecord.builder()
                .from(linx)
                .somaticDrivers(max1(linx.somaticDrivers()))
                .allSomaticStructuralVariants(filterStructuralVariants(linx.allSomaticStructuralVariants(),
                        filteredAllSomaticBreakends, filteredReportableSomaticBreakends))
                .allGermlineStructuralVariants(filterStructuralVariants(linx.allGermlineStructuralVariants(),
                        filteredAllGermlineBreakends, filteredReportableGermlineBreakends))
                .allSomaticFusions(max1(linx.allSomaticFusions()))
                .reportableSomaticFusions(max1(linx.reportableSomaticFusions()))
                .additionalSuspectSomaticFusions(max1(linx.additionalSuspectSomaticFusions()))
                .additionalViableSomaticFusions(max1(linx.additionalViableSomaticFusions()))
                .allSomaticBreakends(filteredAllSomaticBreakends)
                .reportableSomaticBreakends(filteredReportableSomaticBreakends)
                .additionalSuspectSomaticBreakends(max1(linx.additionalSuspectSomaticBreakends()))
                .somaticHomozygousDisruptions(max1(linx.somaticHomozygousDisruptions()))
                .allGermlineStructuralVariants(max1(linx.allGermlineStructuralVariants()))
                .allGermlineBreakends(filteredAllGermlineBreakends)
                .reportableGermlineBreakends(filteredReportableGermlineBreakends)
                .germlineHomozygousDisruptions(max1(linx.germlineHomozygousDisruptions()))
                .build();
    }

    @Nullable
    private static List<LinxSvAnnotation> filterStructuralVariants(@Nullable List<LinxSvAnnotation> allStructuralVariants,
            @Nullable List<LinxBreakend> filteredAllBreakends, @Nullable List<LinxBreakend> filteredReportableBreakends)
    {
        if(allStructuralVariants == null || filteredAllBreakends == null || filteredReportableBreakends == null)
        {
            return null;
        }

        List<Integer> structuralVariantsToRetain = Lists.newArrayList();
        for(LinxBreakend breakend : filteredAllBreakends)
        {
            structuralVariantsToRetain.add(breakend.svId());
        }

        for(LinxBreakend breakend : filteredReportableBreakends)
        {
            structuralVariantsToRetain.add(breakend.svId());
        }

        List<LinxSvAnnotation> filteredStructuralVariants = Lists.newArrayList();
        for(LinxSvAnnotation structuralVariant : allStructuralVariants)
        {
            if(structuralVariantsToRetain.contains(structuralVariant.svId()))
            {
                filteredStructuralVariants.add(structuralVariant);
            }
        }

        return filteredStructuralVariants;
    }

    @Nullable
    private static IsofoxRecord limitIsofoxDataToOne(@Nullable IsofoxRecord isofox)
    {
        if(isofox == null)
        {
            return null;
        }

        return ImmutableIsofoxRecord.builder()
                .from(isofox)
                .allGeneExpressions(max1(isofox.allGeneExpressions()))
                .reportableHighExpression(max1(isofox.reportableHighExpression()))
                .reportableLowExpression(max1(isofox.reportableLowExpression()))
                .allFusions(max1(isofox.allFusions()))
                .reportableNovelKnownFusions(max1(isofox.reportableNovelKnownFusions()))
                .reportableNovelPromiscuousFusions(max1(isofox.reportableNovelPromiscuousFusions()))
                .allNovelSpliceJunctions(max1(isofox.allNovelSpliceJunctions()))
                .reportableSkippedExons(max1(isofox.reportableSkippedExons()))
                .reportableNovelExonsIntrons(max1(isofox.reportableNovelExonsIntrons()))
                .build();
    }

    @Nullable
    private static <T> List<T> max1(@Nullable List<T> elements)
    {
        return elements != null ? elements.subList(0, Math.min(1, elements.size())) : null;
    }
}

