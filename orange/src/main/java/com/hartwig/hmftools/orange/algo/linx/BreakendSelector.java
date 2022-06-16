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

    static final String DOWNSTREAM_ORIENTATION = "Downstream";
    static final String UPSTREAM_ORIENTATION = "Upstream";

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

        KnownFusionData three =
                findByThreeGene(knownFusionCache.getDataByType(KnownFusionType.PROMISCUOUS_3), breakend.gene(), breakend.transcriptId());
        if (three != null) {
            boolean hasReportableFusion = hasReportableThreeFusion(reportableFusions, breakend.gene(), nextExon);
            boolean hasDownstreamOrientation = breakend.geneOrientation().equals(DOWNSTREAM_ORIENTATION);
            boolean isWithinExonRange = isWithinExonRange(three.threeGeneExonRange(), nextExon);
            if (isWithinExonRange && hasDownstreamOrientation && !hasReportableFusion) {
                return true;
            }
        }

        KnownFusionData five =
                findByFiveGene(knownFusionCache.getDataByType(KnownFusionType.PROMISCUOUS_5), breakend.gene(), breakend.transcriptId());
        if (five != null) {
            boolean hasReportableFusion = hasReportableFiveFusion(reportableFusions, breakend.gene(), nextExon);
            boolean hasUpstreamOrientation = breakend.geneOrientation().equals(UPSTREAM_ORIENTATION);
            boolean isWithinExonRange = isWithinExonRange(five.fiveGeneExonRange(), nextExon);

            if (isWithinExonRange && hasUpstreamOrientation && !hasReportableFusion) {
                return true;
            }
        }
        return false;
    }

    private static boolean isWithinExonRange(@NotNull int[] exonRange, int exon) {
        return exon >= exonRange[0] && exon <= exonRange[1];
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
    private static KnownFusionData findByFiveGene(@NotNull List<KnownFusionData> knownFusions, @NotNull String geneToFind,
            @NotNull String transcriptToFind) {
        for (KnownFusionData knownFusion : knownFusions) {
            if (knownFusion.FiveGene.equals(geneToFind) && knownFusion.specificExonsTransName().equals(transcriptToFind)) {
                return knownFusion;
            }
        }
        return null;
    }

    @Nullable
    private static KnownFusionData findByThreeGene(@NotNull List<KnownFusionData> knownFusions, @NotNull String geneToFind,
            @NotNull String transcriptToFind) {
        for (KnownFusionData knownFusion : knownFusions) {
            if (knownFusion.ThreeGene.equals(geneToFind) && knownFusion.specificExonsTransName().equals(transcriptToFind)) {
                return knownFusion;
            }
        }
        return null;
    }
}
