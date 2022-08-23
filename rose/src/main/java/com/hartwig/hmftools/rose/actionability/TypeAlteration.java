package com.hartwig.hmftools.rose.actionability;

import org.jetbrains.annotations.NotNull;

public enum TypeAlteration {
    ACTIVATING_MUTATION,
    ACTIVATING_MUTATION_BRAF_CLASS_I,
    ACTIVATING_MUTATION_BRAF_CLASS_II,
    ACTIVATING_MUTATION_BRAF_CLASS_III,
    ACTIVATING_MUTATION_KRAS_G12C,
    AMPLIFICATION,
    EXTRACELLULAR_DOMAIN_MUTATION,
    FUSION,
    INACTIVATION,
    INTERNAL_DELETION,
    KINASE_DOMAIN_DUPLICATION,
    LOSS,
    POSITIVE,
    RESISTANCE_MUTATION,
    PROMOTER_MUTATION,
    PURITY,
    PURITY_UNRELIABLE,
    NO_ONCOGENIC,
    NO_ACTIONABLE,
    FINDINGS,
    GERMLINE,
    CUPPA,
    CUPPA_INCONCLUSIVE,
    NO_HRD_CAUSE,
    NO_MSI_HRD_PROFILE,
    NOT_BIALLELIC,
    VUS_REMARK,
    UNKNOWN;

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