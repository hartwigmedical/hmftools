package com.hartwig.hmftools.pave.impact;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.VariantData;

public class VariantImpactBuilder
{
    private final GeneDataCache mGeneDataCache;

    public VariantImpactBuilder(final GeneDataCache geneDataCache)
    {
        mGeneDataCache = geneDataCache;
    }

    public VariantImpact createVariantImpact(final VariantData variant)
    {
        // if no impacts were found against any transcripts, then return null and no entry will be written to the VCF
        if(variant.getImpacts().isEmpty())
            return null;

        // find the canonical and worst transcript and their effects
        // favour driver genes over others, otherwise take the gene with the highest impact

        VariantTransImpact worstImpact = null;
        int worstRank = -1;

        VariantTransImpact worstCanonicalImpact = null;
        int worstCanonicalRank = -1;
        String canonicalGeneName = ""; // can be set to the gene name of other reportable transcripts if the canonical has no impact
        String worstGeneName = "";
        String otherReportableTransData = "";

        boolean hasDriverGene = false;

        for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
        {
            final String geneName = entry.getKey();

            boolean isDriverGene = mGeneDataCache.isDriverPanelGene(geneName);

            if(!hasDriverGene && isDriverGene)
            {
                // reset to ignore any prior non-driver gene
                hasDriverGene = true;
                worstImpact = null;
                worstRank = -1;
                worstCanonicalImpact = null;
                worstCanonicalRank = -1;
                worstGeneName = "";
            }
            else if(hasDriverGene && !isDriverGene)
            {
                continue;
            }

            final List<VariantTransImpact> geneImpacts = entry.getValue();

            for(VariantTransImpact transImpact : geneImpacts)
            {
                int rank = transImpact.topRank();

                if(rank > worstRank || (rank == worstRank && transImpact.TransData.IsCanonical))
                {
                    worstRank = rank;
                    worstImpact = transImpact;
                }

                if(transImpact.TransData.IsCanonical && rank > worstCanonicalRank)
                {
                    worstCanonicalRank = rank;
                    worstCanonicalImpact = transImpact;
                    worstGeneName = geneName;
                }
            }

            if(isDriverGene && mGeneDataCache.getOtherReportableTranscripts().containsKey(geneName))
            {
                List<String> otherReportableTrans = mGeneDataCache.getOtherReportableTranscripts().get(geneName);
                StringJoiner sj = new StringJoiner(VAR_IMPACT_OTHER_REPORT_DELIM);

                for(VariantTransImpact transImpact : geneImpacts)
                {
                    if(otherReportableTrans.contains(transImpact.TransData.TransName))
                    {
                        CodingEffect codingEffect = determineCodingEffect(transImpact);

                        sj.add(AltTranscriptReportableInfo.serialise(
                                geneName, transImpact.TransData.TransName, transImpact.hgvsCoding(), transImpact.hgvsProtein(),
                                transImpact.effectsStr(), codingEffect));

                        canonicalGeneName = geneName;
                    }
                }

                otherReportableTransData = sj.toString();
            }
        }

        String canonicalEffect = "";
        String canonicalTranscript = "";
        CodingEffect canonicalCodingEffect = NONE;
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        boolean canonicalSpliceRegion = false;
        CodingEffect worstCodingEffect = NONE;

        if(worstImpact != null)
        {
            worstCodingEffect = determineCodingEffect(worstImpact);
        }

        if(worstCanonicalImpact != null)
        {
            canonicalGeneName = worstGeneName;
            canonicalEffect = worstCanonicalImpact.effectsStr();
            canonicalCodingEffect = determineCodingEffect(worstCanonicalImpact);
            canonicalHgvsCodingImpact = worstCanonicalImpact.hgvsCoding();
            canonicalHgvsProteinImpact = worstCanonicalImpact.hgvsProtein();
            canonicalTranscript = worstCanonicalImpact.TransData.TransName;
            canonicalSpliceRegion = worstCanonicalImpact.inSpliceRegion();
        }

        return new VariantImpact(
                canonicalGeneName, canonicalTranscript, canonicalEffect, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, canonicalSpliceRegion, otherReportableTransData, worstCodingEffect, variant.getImpacts().size());
    }

    private CodingEffect determineCodingEffect(final VariantTransImpact transImpact)
    {
        final List<CodingEffect> simplifiedEffects = transImpact.effects().stream().map(CodingEffect::effect).collect(Collectors.toList());

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT)))
            return NONSENSE_OR_FRAMESHIFT;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE)))
            return SPLICE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE)))
            return MISSENSE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SYNONYMOUS)))
            return SYNONYMOUS;

        return NONE;
    }

}
