package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.SpliceSites.isAcceptorPlusThree;
import static com.hartwig.hmftools.common.variant.SpliceSites.isDonorMinusOne;
import static com.hartwig.hmftools.common.variant.SpliceSites.isDonorPlusFive;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantContextDecorator.getAlt;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class CodingEffectFactory
{
    private final Map<String,List<HmfTranscriptRegion>> mGeneTranscripts; // keyed by transcript ID

    public CodingEffectFactory(final List<HmfTranscriptRegion> transcripts)
    {
        mGeneTranscripts = Maps.newHashMap();

        for(HmfTranscriptRegion transData : transcripts)
        {
            List<HmfTranscriptRegion> geneRegions = mGeneTranscripts.get(transData.geneName());

            if(geneRegions == null)
            {
                geneRegions = Lists.newArrayList();
                mGeneTranscripts.put(transData.geneName(), geneRegions);
            }

            geneRegions.add(transData);
        }
    }

    public boolean hasGene(final String geneName)
    {
        return mGeneTranscripts.containsKey(geneName);
    }

    @NotNull
    public CodingEffect effect(
            @NotNull final VariantContext context, @NotNull final String geneName, final String transName, @NotNull final List<VariantConsequence> consequences)
    {
        final String alt = getAlt(context);
        final VariantType type = VariantType.type(context);

        final List<CodingEffect> simplifiedEffects = consequences.stream().map(CodingEffect::effect).collect(Collectors.toList());

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT)))
            return NONSENSE_OR_FRAMESHIFT;

        List<HmfTranscriptRegion> geneRegions = mGeneTranscripts.get(geneName);

        if(geneRegions != null)
        {
            HmfTranscriptRegion transcript = geneRegions.size() == 1 ?
                    geneRegions.get(0) : geneRegions.stream().filter(x -> x.transName().equals(transName)).findFirst().orElse(null);

            if(transcript != null)
            {
                String acceptorPlusThreeSpliceAlt = transcript.strand() == Strand.FORWARD ? "G" : "C";

                if(consequences.contains(SPLICE_REGION_VARIANT) && (type == VariantType.SNP || type == VariantType.MNP))
                {
                    int position = context.getStart();
                    int end = context.getEnd();

                    while(position <= end)
                    {
                        if(alt.equals(acceptorPlusThreeSpliceAlt) && isAcceptorPlusThree(transcript, position))
                            return SPLICE;

                        if(isDonorMinusOne(transcript, position))
                            return SPLICE;

                        if(isDonorPlusFive(transcript, position))
                            return SPLICE;

                        position++;
                    }
                }
            }
        }

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE)))
            return SPLICE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE)))
            return MISSENSE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SYNONYMOUS)))
            return SYNONYMOUS;

        return NONE;
    }
}
