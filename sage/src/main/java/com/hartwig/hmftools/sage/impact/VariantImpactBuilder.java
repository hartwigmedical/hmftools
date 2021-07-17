package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantImpactBuilder
{
    private final GeneDataCache mGeneDataCache;

    public VariantImpactBuilder(final GeneDataCache geneDataCache)
    {
        mGeneDataCache = geneDataCache;
    }

    public VariantImpact createVariantImapct(final VariantData variant, final VariantContext context)
    {
        String canonicalGene = "";
        String canonicalEffect = "";
        String canonicalTranscript = "";
        CodingEffect canonicalCodingEffect = UNDEFINED;
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        String worstGene = "";
        String worstEffect = "";
        String worstTranscript = "";
        CodingEffect worstCodingEffect = UNDEFINED;

        VariantTransImpact worstImpact = null;
        int topRank = -1;
        Map<String,String> geneIdToNameMap = Maps.newHashMap();

        for(VariantTransImpact transImpact : variant.getImpacts())
        {
            int rank = transImpact.topRank();
            geneIdToNameMap.put(transImpact.TransData.GeneId, findGeneName(transImpact.TransData));

            if(rank > topRank || (rank == topRank && transImpact.TransData.IsCanonical))
            {
                topRank = rank;
                worstImpact = transImpact;
            }
        }

        if(worstImpact != null)
        {
            worstGene = geneIdToNameMap.get(worstImpact.TransData.GeneId);
            worstEffect = VariantConsequence.consequenceString(worstImpact.consequences());
            worstCodingEffect = determineCodingEffect(variant, worstImpact, context);
            worstTranscript = worstImpact.TransData.TransName;
        }

        List<VariantTransImpact> canonicalImpacts = variant.getImpacts().stream()
                .filter(x -> x.TransData.IsCanonical).collect(Collectors.toList());

        List<VariantTransImpact> canonicalDriverImpacts = canonicalImpacts.stream()
                .filter(x -> mGeneDataCache.getDriverPanelGenes().contains(findGeneName(x.TransData)))
                .collect(Collectors.toList());

        if(!canonicalImpacts.isEmpty())
        {
            VariantTransImpact canonicalImpact = !canonicalDriverImpacts.isEmpty() ?
                    canonicalDriverImpacts.get(0) : canonicalImpacts.get(0);

            canonicalGene = geneIdToNameMap.get(canonicalImpact.TransData.GeneId);
            canonicalEffect = VariantConsequence.consequenceString(canonicalImpact.consequences());
            canonicalCodingEffect = determineCodingEffect(variant, canonicalImpact, context);
            canonicalHgvsCodingImpact = canonicalImpact.hgvsCodingChange();
            canonicalHgvsProteinImpact = canonicalImpact.hgvsProteinChange();
            canonicalTranscript = canonicalImpact.TransData.TransName;
        }

        return new VariantImpact(
                geneIdToNameMap.size(), canonicalGene, canonicalEffect, canonicalTranscript, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, worstGene, worstEffect, worstTranscript, worstCodingEffect);
    }

    private CodingEffect determineCodingEffect(final VariantData variant, final VariantTransImpact transImpact, final VariantContext context)
    {
        final List<CodingEffect> simplifiedEffects = transImpact.consequences().stream().map(CodingEffect::effect).collect(Collectors.toList());

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT)))
        {
            boolean phasedInframeIndel = context.isIndel() && context.getAttributeAsInt(SageMetaData.PHASED_INFRAME_INDEL, 0) > 0;

            if(phasedInframeIndel)
                return MISSENSE;

            return NONSENSE_OR_FRAMESHIFT;
        }

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE)))
            return SPLICE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE)))
            return MISSENSE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SYNONYMOUS)))
            return SYNONYMOUS;

        return NONE;
    }

    private String findGeneName(final TranscriptData transData)
    {
        return mGeneDataCache.getEnsemblCache().getGeneDataByName(transData.GeneId).GeneName;
    }

}
