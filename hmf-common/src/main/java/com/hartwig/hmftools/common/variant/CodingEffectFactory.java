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

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class CodingEffectFactory {

    private final Map<String, HmfTranscriptRegion> transcripts;

    public CodingEffectFactory(final List<HmfTranscriptRegion> transcripts) {
        this.transcripts = transcripts.stream().collect(Collectors.toMap(TranscriptRegion::gene, x -> x));
    }

    @NotNull
    public CodingEffect effect(
            @NotNull final VariantContext context, @NotNull final String gene, @NotNull final List<VariantConsequence> consequences) {
        final VariantContextDecorator variant = new VariantContextDecorator(context);
        final String alt = variant.alt();

        final List<CodingEffect> simplifiedEffects = consequences.stream().map(CodingEffect::effect).collect(Collectors.toList());

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT))) {
            return NONSENSE_OR_FRAMESHIFT;
        }

        HmfTranscriptRegion transcript = transcripts.get(gene);
        if (transcript != null) {

            if (variant.type() == VariantType.SNP || variant.type() == VariantType.MNP) {
                int position = variant.context().getStart();
                int end = variant.context().getEnd();
                while (position <= end) {
                    String acceptorPlusThreeSpliceAlt = transcript.strand() == Strand.FORWARD ? "G" : "C";
                    if (consequences.contains(SPLICE_REGION_VARIANT) && alt.equals(acceptorPlusThreeSpliceAlt) && isAcceptorPlusThree(transcript, position)) {
                        return SPLICE;
                    }

                    if (consequences.contains(SPLICE_REGION_VARIANT) && isDonorMinusOne(transcript, position)) {
                        return SPLICE;
                    }

                    if (consequences.contains(SPLICE_REGION_VARIANT) && isDonorPlusFive(transcript, position)) {
                        return SPLICE;
                    }

                    position++;
                }
            }

        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE))) {
            return SPLICE;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE))) {
            return MISSENSE;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(SYNONYMOUS))) {
            return SYNONYMOUS;
        }

        return NONE;
    }

}
