package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;

public enum CodingEffect {
    NONSENSE_OR_FRAMESHIFT,
    SPLICE,
    MISSENSE,
    SYNONYMOUS,
    NONE,
    UNDEFINED;

    @NotNull
    public static CodingEffect effect(@NotNull final String gene, @NotNull final List<VariantConsequence> consequences) {
        final List<CodingEffect> simplifiedEffects = consequences.stream().map(CodingEffect::effect).collect(Collectors.toList());

        if (gene.equals("TP53")) {
            // TP53 has some known pathogenic variants in splice regions on the exon-side
            if (consequences.contains(SYNONYMOUS_VARIANT) && consequences.contains(SPLICE_REGION_VARIANT)) {
                return SPLICE;
            }
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT))) {
            return NONSENSE_OR_FRAMESHIFT;
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

    @NotNull
    public static CodingEffect effect(@NotNull final HmfTranscriptRegion transcript, long position, String alt,
            @NotNull final List<VariantConsequence> consequences) {
        final List<CodingEffect> simplifiedEffects = consequences.stream().map(CodingEffect::effect).collect(Collectors.toList());

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT))) {
            return NONSENSE_OR_FRAMESHIFT;
        }

        if (consequences.contains(SPLICE_REGION_VARIANT) && alt.equals("G") && transcript.isAcceptorPlusThree(position)) {
            return SPLICE;
        }

        if (consequences.contains(SPLICE_REGION_VARIANT) && transcript.isDonorMinusOne(position)) {
            return SPLICE;
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

    @NotNull
    private static CodingEffect effect(@NotNull final VariantConsequence consequence) {
        switch (consequence) {
            case FRAMESHIFT_VARIANT:
            case STOP_GAINED:
            case STOP_LOST:
            case START_LOST:
            case INITIATOR_CODON_VARIANT:
                return NONSENSE_OR_FRAMESHIFT;
            case MISSENSE_VARIANT:
            case PROTEIN_PROTEIN_CONTACT:
            case STRUCTURAL_INTERACTION_VARIANT:
            case INFRAME_DELETION:
            case INFRAME_INSERTION:
                return MISSENSE;
            case SPLICE_ACCEPTOR_VARIANT:
            case SPLICE_DONOR_VARIANT:
                return SPLICE;
            case SYNONYMOUS_VARIANT:
                return SYNONYMOUS;
        }

        return NONE;
    }
}
