package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.common.rna.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.rna.KnownFusionType.hasPromiscousGene;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;

import org.jetbrains.annotations.Nullable;

public final class RnaFusionSelector
{
    private static final Set<StructuralVariantType> ALWAYS_VALID_FOR_PROMISCUOUS = Sets.newHashSet();

    static
    {
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.BND);
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.INV);
        ALWAYS_VALID_FOR_PROMISCUOUS.add(StructuralVariantType.INS);
    }

    private static final String RNA_FUSION_NAME_DELIMITER = "_";

    @Nullable
    public static String geneUp(final RnaFusion rnaFusion)
    {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split > 0 ? rnaFusion.name().substring(0, split) : null;
    }

    @Nullable
    public static String geneDown(final RnaFusion rnaFusion)
    {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split >= 0 && split < rnaFusion.name().length() - 1 ? rnaFusion.name().substring(split + 1) : null;
    }

    public static boolean hasFusion(final List<LinxFusion> linxFusions, final String geneUp, final String geneDown)
    {
        for(LinxFusion linxFusion : linxFusions)
        {
            if(linxFusion.geneStart().equals(geneUp) && linxFusion.geneEnd().equals(geneDown))
            {
                return true;
            }
        }
        return false;
    }

    public static List<RnaFusion> selectNovelKnownFusions(final List<RnaFusion> rnaFusions, final List<LinxFusion> linxFusions)
    {
        List<RnaFusion> result = Lists.newArrayList();

        for(RnaFusion rnaFusion : rnaFusions)
        {
            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if(geneUp != null && geneDown != null)
            {
                if(rnaFusion.knownType() == KNOWN_PAIR && !hasFusion(linxFusions, geneUp, geneDown))
                {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }

    public static List<RnaFusion> selectNovelPromiscuousFusions(final List<RnaFusion> rnaFusions, final List<LinxFusion> linxFusions)
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
                boolean isPromiscuous = hasPromiscousGene(rnaFusion.knownType());
                boolean isKnown = rnaFusion.knownType() == KNOWN_PAIR;
                if(isPromiscuous && !isKnown && !hasFusion(linxFusions, geneUp, geneDown))
                {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }
}
