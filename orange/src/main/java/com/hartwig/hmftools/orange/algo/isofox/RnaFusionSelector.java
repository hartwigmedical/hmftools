package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneDown;
import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneUp;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.finding.Fusion;
import com.hartwig.hmftools.orange.algo.linx.DnaFusionEvaluator;

import org.jetbrains.annotations.NotNull;

final class RnaFusionSelector
{
    private static final Set<StructuralVariantType> ALWAYS_VALID_FOR_PROMISCUOUS = Sets.newHashSet();

    static
    {
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.BND);
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.INV);
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.INS);
    }

    @NotNull
    public static List<RnaFusion> selectNovelKnownFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<Fusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache)
    {
        List<RnaFusion> result = Lists.newArrayList();

        for(RnaFusion rnaFusion : rnaFusions)
        {
            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if(geneUp != null && geneDown != null)
            {
                if(knownFusionCache.hasKnownFusion(geneUp, geneDown) && !DnaFusionEvaluator.hasFusion(linxFusions, geneUp, geneDown))
                {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }

    @NotNull
    public static List<RnaFusion> selectNovelPromiscuousFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<Fusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache)
    {
        List<RnaFusion> result = Lists.newArrayList();

        for(RnaFusion rnaFusion : rnaFusions)
        {
            boolean isTypeMatch = ALWAYS_VALID_FOR_PROMISCUOUS.contains(rnaFusion.svType());
            boolean hasSufficientDistance = Math.abs(rnaFusion.positionUp() - rnaFusion.positionDown()) > 1E6;

            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if(geneUp != null && geneDown != null && (isTypeMatch || hasSufficientDistance))
            {
                boolean isPromiscuous =
                        knownFusionCache.hasPromiscuousFiveGene(geneUp) || knownFusionCache.hasPromiscuousThreeGene(geneDown);
                boolean isKnown = knownFusionCache.hasKnownFusion(geneUp, geneDown);
                if(isPromiscuous && !isKnown && !DnaFusionEvaluator.hasFusion(linxFusions, geneUp, geneDown))
                {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }
}
