package com.hartwig.hmftools.summon.actionability;

import org.jetbrains.annotations.NotNull;

public enum TypeAlteration {
    ACTIVATING_MUTATION,
    AMPLIFICATION,
    EXTRACELLULAR_DOMAIN_MUTATION,
    FUSION,
    INACTIVATION,
    INTERNAL_DELETION,
    KINASE_DOMAIN_DUPLICATION,
    LOSS,
    POSITIVE,
    RESISTANCE_MUTATION,
    PURITY,
    NO_ONCOGENIC,
    NO_ACTIONABLE,
    FINDINGS,
    GERMLINE,
    CUPPA,
    CUPPA_INCONCLUSIVE,
    NO_HRD_CAUSE,
    NO_MSI_HRD_PROFILE,
    NOT_BIALLELIC;

    @NotNull
    static TypeAlteration toType(@NotNull String typeInput) {
        for (TypeAlteration type : TypeAlteration.values()) {
            if (typeInput.equals(type.toString())) {
                return type;
            }
        }

        throw new IllegalStateException("Cannot resolve type: '{}' " + typeInput);
    }
}