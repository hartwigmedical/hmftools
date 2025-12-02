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

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public final class ReportLimiter
{
    public static OrangeRecord limitAllListsToMaxOne(final OrangeRecord report)
    {
        PurpleRecord purple = limitPurpleDataToOne(report.purple());
        LinxRecord linx = limitLinxDataToOne(report.linx());

        return ImmutableOrangeRecord.builder()
                .from(report)
                .germlineMVLHPerGene(limitGermlineMVLHToOne(report.germlineMVLHPerGene()))
                .purple(purple)
                .linx(linx)
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

    private static PurpleRecord limitPurpleDataToOne(final PurpleRecord purple)
    {
        return ImmutablePurpleRecord.builder()
                .from(purple)
                .somaticDrivers(max1(purple.somaticDrivers()))
                .germlineDrivers(max1(purple.germlineDrivers()))
                .otherSomaticVariants(max1(purple.otherSomaticVariants()))
                .driverSomaticVariants(max1(purple.driverSomaticVariants()))
                .otherGermlineVariants(max1(purple.otherGermlineVariants()))
                .driverGermlineVariants(max1(purple.driverGermlineVariants()))
                .somaticCopyNumbers(max1(purple.somaticCopyNumbers()))
                .somaticGeneCopyNumbers(max1(purple.somaticGeneCopyNumbers()))
                .driverSomaticGainsDels(max1(purple.driverSomaticGainsDels()))
                .otherGermlineDeletions(max1(purple.otherGermlineDeletions()))
                .driverGermlineDeletions(max1(purple.driverGermlineDeletions()))
                .allGermlineLossOfHeterozygosities(max1(purple.allGermlineLossOfHeterozygosities()))
                .driverGermlineLossOfHeterozygosities(max1(purple.driverGermlineLossOfHeterozygosities()))
                .build();
    }

    private static LinxRecord limitLinxDataToOne(final LinxRecord linx)
    {
        List<LinxBreakend> filteredAllSomaticBreakends = max1(linx.otherSomaticBreakends());
        List<LinxBreakend> filteredReportableSomaticBreakends = max1(linx.driverSomaticBreakends());
        List<LinxBreakend> filteredAllGermlineBreakends = max1(linx.otherGermlineBreakends());
        List<LinxBreakend> filteredReportableGermlineBreakends = max1(linx.driverGermlineBreakends());

        return ImmutableLinxRecord.builder()
                .from(linx)
                .somaticDrivers(max1(linx.somaticDrivers()))
                .allSomaticStructuralVariants(filterStructuralVariants(linx.allSomaticStructuralVariants(),
                        filteredAllSomaticBreakends, filteredReportableSomaticBreakends))
                .allGermlineStructuralVariants(filterStructuralVariants(linx.allGermlineStructuralVariants(),
                        filteredAllGermlineBreakends, filteredReportableGermlineBreakends))
                .allSomaticFusions(max1(linx.allSomaticFusions()))
                .otherSomaticBreakends(filteredAllSomaticBreakends)
                .driverSomaticBreakends(filteredReportableSomaticBreakends)
                .somaticHomozygousDisruptions(max1(linx.somaticHomozygousDisruptions()))
                .allGermlineStructuralVariants(max1(linx.allGermlineStructuralVariants()))
                .otherGermlineBreakends(filteredAllGermlineBreakends)
                .driverGermlineBreakends(filteredReportableGermlineBreakends)
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

