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

    private static PurpleRecord limitPurpleDataToOne(final PurpleRecord purple)
    {
        return ImmutablePurpleRecord.builder()
                .from(purple)
                .somaticDrivers(max1(purple.somaticDrivers()))
                .germlineDrivers(max1(purple.germlineDrivers()))
                .somaticVariants(max1(purple.somaticVariants()))
                .germlineVariants(max1(purple.germlineVariants()))
                .somaticCopyNumbers(max1(purple.somaticCopyNumbers()))
                .somaticGeneCopyNumbers(max1(purple.somaticGeneCopyNumbers()))
                .somaticGainsDels(max1(purple.somaticGainsDels()))
                .germlineGainsDels(max1(purple.germlineGainsDels()))
                .build();
    }

    private static LinxRecord limitLinxDataToOne(final LinxRecord linx)
    {
        List<LinxBreakend> filteredSomaticBreakends = max1(linx.somaticBreakends());
        List<LinxBreakend> filteredGermlineBreakends = max1(linx.germlineBreakends());

        return ImmutableLinxRecord.builder()
                .from(linx)
                .somaticDrivers(max1(linx.somaticDrivers()))
                .somaticStructuralVariants(filterStructuralVariants(linx.somaticStructuralVariants(), filteredSomaticBreakends))
                .germlineStructuralVariants(filterStructuralVariants(linx.germlineStructuralVariants(), filteredGermlineBreakends))
                .fusions(max1(linx.fusions()))
                .somaticBreakends(filteredSomaticBreakends)
                .somaticHomozygousDisruptions(max1(linx.somaticHomozygousDisruptions()))
                .germlineBreakends(filteredGermlineBreakends)
                .germlineHomozygousDisruptions(max1(linx.germlineHomozygousDisruptions()))
                .build();
    }

    @Nullable
    private static List<LinxSvAnnotation> filterStructuralVariants(
            @Nullable List<LinxSvAnnotation> structuralVariants, @Nullable List<LinxBreakend> filteredBreakends)
    {
        if(structuralVariants == null || filteredBreakends == null)
        {
            return null;
        }

        List<Integer> structuralVariantsToRetain = Lists.newArrayList();
        for(LinxBreakend breakend : filteredBreakends)
        {
            structuralVariantsToRetain.add(breakend.svId());
        }

        List<LinxSvAnnotation> filteredStructuralVariants = Lists.newArrayList();
        for(LinxSvAnnotation structuralVariant : structuralVariants)
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

