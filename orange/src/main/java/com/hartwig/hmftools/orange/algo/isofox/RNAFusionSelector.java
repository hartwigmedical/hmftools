package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.orange.algo.linx.DNAFusionEvaluator;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class RNAFusionSelector {

    private static final String RNA_FUSION_NAME_DELIMITER = "_";
    private static final Set<StructuralVariantType> ALWAYS_VALID_FOR_PROMISCUOUS = Sets.newHashSet();

    static {
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.BND);
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.INV);
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.INS);
    }

    private RNAFusionSelector() {
    }

    @NotNull
    public static List<RnaFusion> selectNovelKnownFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<LinxFusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache) {
        List<RnaFusion> result = Lists.newArrayList();

        for (RnaFusion rnaFusion : rnaFusions) {
            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if (geneUp != null && geneDown != null) {
                if (knownFusionCache.hasKnownFusion(geneUp, geneDown) && !DNAFusionEvaluator.hasFusion(linxFusions, geneUp, geneDown)) {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }

    @NotNull
    public static List<RnaFusion> selectNovelPromiscuousFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<LinxFusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache) {
        List<RnaFusion> result = Lists.newArrayList();

        for (RnaFusion rnaFusion : rnaFusions) {
            boolean isTypeMatch = ALWAYS_VALID_FOR_PROMISCUOUS.contains(rnaFusion.svType());
            boolean hasSufficientDistance = Math.abs(rnaFusion.positionUp() - rnaFusion.positionDown()) > 1E6;

            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if (geneUp != null && geneDown != null && (isTypeMatch || hasSufficientDistance)) {
                boolean isPromiscuous =
                        knownFusionCache.hasPromiscuousFiveGene(geneUp) || knownFusionCache.hasPromiscuousThreeGene(geneDown);
                boolean isKnown = knownFusionCache.hasKnownFusion(geneUp, geneDown);
                if (isPromiscuous && !isKnown && !DNAFusionEvaluator.hasFusion(linxFusions, geneUp, geneDown)) {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }

    @VisibleForTesting
    @Nullable
    static String geneUp(@NotNull RnaFusion rnaFusion) {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split > 0 ? rnaFusion.name().substring(0, split) : null;
    }

    @VisibleForTesting
    @Nullable
    static String geneDown(@NotNull RnaFusion rnaFusion) {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split >= 0 && split < rnaFusion.name().length() - 1 ? rnaFusion.name().substring(split + 1) : null;
    }
}
