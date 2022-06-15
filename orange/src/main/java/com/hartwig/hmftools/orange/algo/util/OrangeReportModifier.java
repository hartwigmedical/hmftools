package com.hartwig.hmftools.orange.algo.util;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.isofox.ImmutableIsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.ImmutableLinxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;

import org.jetbrains.annotations.NotNull;

public final class OrangeReportModifier {

    private OrangeReportModifier() {
    }

    @NotNull
    public static OrangeReport limitAllListsToMaxOne(@NotNull OrangeReport report) {
        return ImmutableOrangeReport.builder()
                .from(report)
                .germlineMVLHPerGene(limitGermlineMVLHToOne(report.germlineMVLHPerGene()))
                .purple(limitPurpleDataToOne(report.purple()))
                .linx(limitLinxDataToOne(report.linx()))
                .isofox(limitIsofoxDataToOne(report.isofox()))
                .protect(limitProtectDataToOne(report.protect()))
                .build();
    }

    @NotNull
    private static Map<String, Double> limitGermlineMVLHToOne(@NotNull Map<String, Double> germlineMVLHPerGene) {
        Map<String, Double> filtered = Maps.newHashMap();
        if (!germlineMVLHPerGene.isEmpty()) {
            String firstKey = germlineMVLHPerGene.keySet().iterator().next();
            filtered.put(firstKey, germlineMVLHPerGene.get(firstKey));
        }
        return filtered;
    }

    @NotNull
    private static PurpleInterpretedData limitPurpleDataToOne(@NotNull PurpleInterpretedData purple) {
        return ImmutablePurpleInterpretedData.builder()
                .from(purple)
                .allSomaticVariants(max1(purple.allSomaticVariants()))
                .reportableSomaticVariants(max1(purple.reportableSomaticVariants()))
                .additionalSuspectSomaticVariants(max1(purple.additionalSuspectSomaticVariants()))
                .allGermlineVariants(max1(purple.allGermlineVariants()))
                .reportableGermlineVariants(max1(purple.reportableGermlineVariants()))
                .additionalSuspectGermlineVariants(max1(purple.additionalSuspectGermlineVariants()))
                .allSomaticGeneCopyNumbers(max1(purple.allSomaticGeneCopyNumbers()))
                .suspectGeneCopyNumbersWithLOH(max1(purple.suspectGeneCopyNumbersWithLOH()))
                .allSomaticGainsLosses(max1(purple.allSomaticGainsLosses()))
                .reportableSomaticGainsLosses(max1(purple.reportableSomaticGainsLosses()))
                .nearReportableSomaticGains(max1(purple.nearReportableSomaticGains()))
                .additionalSuspectSomaticGainsLosses(max1(purple.additionalSuspectSomaticGainsLosses()))
                .allGermlineDeletions(max1(purple.allGermlineDeletions()))
                .reportableGermlineDeletions(max1(purple.reportableGermlineDeletions()))
                .copyNumberPerChromosome(max1(purple.copyNumberPerChromosome()))
                .build();
    }

    @NotNull
    private static LinxInterpretedData limitLinxDataToOne(@NotNull LinxInterpretedData linx) {
        return ImmutableLinxInterpretedData.builder()
                .from(linx)
                .allFusions(max1(linx.allFusions()))
                .reportableFusions(max1(linx.reportableFusions()))
                .additionalSuspectFusions(max1(linx.additionalSuspectFusions()))
                .allBreakends(max1(linx.allBreakends()))
                .reportableGeneDisruptions(max1(linx.reportableGeneDisruptions()))
                .additionalSuspectBreakends(max1(linx.additionalSuspectBreakends()))
                .homozygousDisruptions(max1(linx.homozygousDisruptions()))
                .allGermlineDisruptions(max1(linx.allGermlineDisruptions()))
                .reportableGermlineDisruptions(max1(linx.reportableGermlineDisruptions()))
                .build();
    }

    @NotNull
    private static IsofoxInterpretedData limitIsofoxDataToOne(@NotNull IsofoxInterpretedData isofox) {
        return ImmutableIsofoxInterpretedData.builder()
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

    @NotNull
    private static List<ProtectEvidence> limitProtectDataToOne(@NotNull List<ProtectEvidence> protect) {
        return max1(protect);
    }

    @NotNull
    private static <T> List<T> max1(@NotNull List<T> elements) {
        return elements.subList(0, Math.min(1, elements.size()));
    }
}

