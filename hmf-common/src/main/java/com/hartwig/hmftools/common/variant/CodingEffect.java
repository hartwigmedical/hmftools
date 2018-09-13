package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.stream.Collectors;

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
        final List<CodingEffect> simplifiedEffects = consequences.stream().map(x -> effect(gene, x)).collect(Collectors.toList());
        if (simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT))) {
            return NONSENSE_OR_FRAMESHIFT;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE))) {
            return MISSENSE;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE))) {
            return SPLICE;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(SYNONYMOUS))) {
            return SYNONYMOUS;
        }

        return NONE;
    }

    @NotNull
    private static CodingEffect effect(@NotNull final String gene, @NotNull final VariantConsequence consequence) {
        // KODU: Below exception exists because TP53 has some known pathogenic variants in splice regions.
        if (gene.equals("TP53") && consequence.equals(VariantConsequence.SPLICE_REGION_VARIANT)) {
            return SPLICE;
        }

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
