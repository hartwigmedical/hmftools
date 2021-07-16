package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.isAcceptorPlusThree;
import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.isDonorMinusOne;
import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.isDonorPlusFive;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantContextDecorator.getAlt;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class CodingEffectFactory
{
    private final Map<String, HmfTranscriptRegion> mTranscripts;

    public CodingEffectFactory(final List<HmfTranscriptRegion> transcripts)
    {
        mTranscripts = transcripts.stream().collect(Collectors.toMap(TranscriptRegion::gene, x -> x));
    }

    @NotNull
    public CodingEffect effect(
            @NotNull final VariantContext context, @NotNull final String gene, @NotNull final List<VariantConsequence> consequences)
    {
        final String alt = getAlt(context);
        final VariantType type = VariantType.type(context);

        final List<CodingEffect> simplifiedEffects = consequences.stream().map(CodingEffect::effect).collect(Collectors.toList());

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT)))
            return NONSENSE_OR_FRAMESHIFT;

        HmfTranscriptRegion transcript = mTranscripts.get(gene);
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

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE)))
            return SPLICE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE)))
            return MISSENSE;

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SYNONYMOUS)))
            return SYNONYMOUS;

        return NONE;
    }
}
