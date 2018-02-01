package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public enum SimplifiedEffect {
    NONSENSE,
    SPLICE,
    MISSENSE,
    INDEL_FRAMESHIFT,
    INDEL_OTHER,
    NONE;

    @NotNull
    public static SimplifiedEffect effect(@NotNull final VariantType type, @NotNull final List<VariantConsequence> consequences) {
        switch (type) {
            case INDEL:
                return indelEffect(consequences);
            case SNP:
            case MNP:
                return snpEffect(consequences);
        }

        return NONE;
    }

    @NotNull
    public static SimplifiedEffect indelEffect(@NotNull final List<VariantConsequence> consequences) {

        final List<SimplifiedEffect> simplifiedEffects =
                consequences.stream().map(SimplifiedEffect::indelEffect).collect(Collectors.toList());
        if (simplifiedEffects.stream().anyMatch(x -> x.equals(INDEL_FRAMESHIFT))) {
            return INDEL_FRAMESHIFT;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(INDEL_OTHER))) {
            return INDEL_OTHER;
        }

        return NONE;
    }

    @NotNull
    public static SimplifiedEffect snpEffect(@NotNull final List<VariantConsequence> consequences) {

        final List<SimplifiedEffect> simplifiedEffects =
                consequences.stream().map(SimplifiedEffect::snpEffect).collect(Collectors.toList());
        if (simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE))) {
            return NONSENSE;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE))) {
            return SPLICE;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE))) {
            return MISSENSE;
        }

        return NONE;
    }

    @NotNull
    private static SimplifiedEffect snpEffect(@NotNull final VariantConsequence consequence) {

        if (!consequence.isActionable()) {
            return NONE;
        }

        switch (consequence) {
            case INTRON_VARIANT:
            case UTR_VARIANT:
            case SEQUENCE_FEATURE:
            case NON_CODING_TRANSCRIPT_VARIANT:
            case SYNONYMOUS_VARIANT:
                return NONE;

            case STOP_GAINED:
                return NONSENSE;
            case SPLICE_ACCEPTOR_VARIANT:
            case SPLICE_DONOR_VARIANT:
            case SPLICE_REGION_VARIANT:
                return SPLICE;
        }

        return MISSENSE;
    }

    @NotNull
    private static SimplifiedEffect indelEffect(@NotNull final VariantConsequence consequence) {

        if (!consequence.isActionable()) {
            return NONE;
        }

        switch (consequence) {
            case STOP_GAINED:
            case FRAMESHIFT_VARIANT:
                return INDEL_FRAMESHIFT;
            case PROTEIN_PROTEIN_CONTACT:
            case STRUCTURAL_INTERACTION_VARIANT:
            case INFRAME_DELETION:
            case INFRAME_INSERTION:
                return INDEL_OTHER;
        }

        return MISSENSE;
    }

}
