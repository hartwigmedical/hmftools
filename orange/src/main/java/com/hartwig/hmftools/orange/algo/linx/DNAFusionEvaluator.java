package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;

public final class DNAFusionEvaluator {

    private DNAFusionEvaluator() {
    }

    public static boolean hasFusion(@NotNull List<LinxFusion> linxFusions, @NotNull String geneUp, @NotNull String geneDown) {
        for (LinxFusion linxFusion : linxFusions) {
            if (linxFusion.geneStart().equals(geneUp) && linxFusion.geneEnd().equals(geneDown)) {
                return true;
            }
        }
        return false;
    }
}
