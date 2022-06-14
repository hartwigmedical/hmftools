package com.hartwig.hmftools.orange.util;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.isofox.ImmutableIsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;

import org.jetbrains.annotations.NotNull;

public final class OrangeReportModifier {

    private OrangeReportModifier() {
    }

    @NotNull
    public static OrangeReport limitAllListsToMaxOne(@NotNull OrangeReport report) {
        return ImmutableOrangeReport.builder().from(report)
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
    private static PurpleData limitPurpleDataToOne(@NotNull PurpleData purple) {
        return ImmutablePurpleData.builder().from(purple)
                .reportableSomaticVariants(max1(purple.reportableSomaticVariants()))
                .unreportedSomaticVariants(max1(purple.unreportedSomaticVariants()))
                .reportableGermlineVariants(max1(purple.reportableGermlineVariants()))
                .unreportedGermlineVariants(max1(purple.unreportedGermlineVariants()))
                .reportableGermlineDeletions(max1(purple.reportableGermlineDeletions()))
                .unreportedGermlineDeletions(max1(purple.unreportedGermlineDeletions()))
                .allGeneCopyNumbers(max1(purple.allGeneCopyNumbers()))
                .copyNumberPerChromosome(max1(purple.copyNumberPerChromosome()))
                .build();
    }

    @NotNull
    private static LinxData limitLinxDataToOne(@NotNull LinxData linx) {
        return ImmutableLinxData.builder().reportableFusions(max1(linx.reportableFusions()))
                .unreportedFusions(max1(linx.unreportedFusions()))
                .geneDisruptions(max1(linx.geneDisruptions()))
                .homozygousDisruptions(max1(linx.homozygousDisruptions()))
                .drivers(max1(linx.drivers()))
                .reportableGermlineDisruptions(max1(linx.reportableGermlineDisruptions()))
                .unreportedGermlineDisruptions(max1(linx.unreportedGermlineDisruptions()))
                .build();
    }

    @NotNull
    private static IsofoxInterpretedData limitIsofoxDataToOne(@NotNull IsofoxInterpretedData isofox) {
        return ImmutableIsofoxInterpretedData.builder().from(isofox)
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

