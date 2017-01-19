package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class VariantAnalysis {

    @NotNull
    private final List<SomaticVariant> allVariants;
    @NotNull
    private final List<SomaticVariant> allPassedVariants;
    @NotNull
    private final List<SomaticVariant> consensusPassedVariants;
    @NotNull
    private final List<SomaticVariant> missenseVariants;
    @NotNull
    private final List<SomaticVariant> consequencePassedVariants;

    VariantAnalysis(@NotNull final List<SomaticVariant> allVariants,
            @NotNull final List<SomaticVariant> allPassedVariants,
            @NotNull final List<SomaticVariant> consensusPassedVariants,
            @NotNull final List<SomaticVariant> missenseVariants,
            @NotNull final List<SomaticVariant> consequencePassedVariants) {
        this.allVariants = allVariants;
        this.allPassedVariants = allPassedVariants;
        this.consensusPassedVariants = consensusPassedVariants;
        this.missenseVariants = missenseVariants;
        this.consequencePassedVariants = consequencePassedVariants;
    }

    @NotNull
    public List<SomaticVariant> allVariants() {
        return allVariants;
    }

    @NotNull
    public List<SomaticVariant> allPassedVariants() {
        return allPassedVariants;
    }

    @NotNull
    public List<SomaticVariant> consensusPassedVariants() {
        return consensusPassedVariants;
    }

    @NotNull
    public List<SomaticVariant> missenseVariants() {
        return missenseVariants;
    }

    @NotNull
    public List<SomaticVariant> consequencePassedVariants() {
        return consequencePassedVariants;
    }
}
