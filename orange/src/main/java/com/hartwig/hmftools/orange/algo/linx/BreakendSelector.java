package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class BreakendSelector {

    private BreakendSelector() {
    }

    @NotNull
    public static List<LinxBreakend> selectInterestingUnreportedBreakends(@NotNull List<LinxBreakend> allBreakends,
            @NotNull List<LinxFusion> reportableFusions, @NotNull KnownFusionCache knownFusionCache) {
        List<LinxBreakend> interestingUnreportedBreakends = Lists.newArrayList();
        for (LinxBreakend breakend : allBreakends) {
            if (!breakend.reportedDisruption() && breakend.disruptive()) {
                if (isUnreportedBreakInPromiscuousExonRange(knownFusionCache, reportableFusions, breakend)) {
                    interestingUnreportedBreakends.add(breakend);
                }
            }
        }
        return interestingUnreportedBreakends;
    }

    private static boolean isUnreportedBreakInPromiscuousExonRange(@NotNull KnownFusionCache knownFusionCache,
            @NotNull List<LinxFusion> reportableFusions, @NotNull LinxBreakend breakend) {
        int nextExon = breakend.nextSpliceExonRank();

        KnownFusionData three = findByThreeGene(knownFusionCache.getDataByType(KnownFusionType.PROMISCUOUS_3), breakend.gene());
        if (three != null) {
            boolean hasReportedFusion = hasReportableThreeFusion(reportableFusions, breakend.gene(), nextExon);
            int[] exonRange = three.threeGeneExonRange();
            if (nextExon >= exonRange[0] && nextExon <= exonRange[1] && !hasReportedFusion) {
                return true;
            }
        }

        KnownFusionData five = findByThreeGene(knownFusionCache.getDataByType(KnownFusionType.PROMISCUOUS_5), breakend.gene());
        if (five != null) {
            boolean hasReportedFusion = hasReportableFiveFusion(reportableFusions, breakend.gene(), nextExon);
            int[] exonRange = three.fiveGeneExonRange();
            if (nextExon >= exonRange[0] && nextExon <= exonRange[1] && !hasReportedFusion) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasReportableFiveFusion(@NotNull List<LinxFusion> reportableFusions, @NotNull String gene, int exon) {
        for (LinxFusion fusion : reportableFusions) {
            if (fusion.geneStart().equals(gene) && fusion.fusedExonUp() == exon) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasReportableThreeFusion(@NotNull List<LinxFusion> reportableFusions, @NotNull String gene, int exon) {
        for (LinxFusion fusion : reportableFusions) {
            if (fusion.geneEnd().equals(gene) && fusion.fusedExonDown() == exon) {
                return true;
            }
        }
        return false;
    }

    @Nullable
    private static KnownFusionData findByFiveGene(@NotNull List<KnownFusionData> knownFusions, @NotNull String geneToFind) {
        for (KnownFusionData knownFusion : knownFusions) {
            if (knownFusion.FiveGene.equals(geneToFind)) {
                return knownFusion;
            }
        }
        return null;
    }

    @Nullable
    private static KnownFusionData findByThreeGene(@NotNull List<KnownFusionData> knownFusions, @NotNull String geneToFind) {
        for (KnownFusionData knownFusion : knownFusions) {
            if (knownFusion.ThreeGene.equals(geneToFind)) {
                return knownFusion;
            }
        }
        return null;
    }
}
