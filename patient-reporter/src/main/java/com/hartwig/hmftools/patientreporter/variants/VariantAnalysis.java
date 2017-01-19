package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class VariantAnalysis {

    @NotNull
    private final List<SomaticVariant> passedVariants;
    @NotNull
    private final List<SomaticVariant> consensusPassedVariants;
    @NotNull
    private final List<SomaticVariant> missenseVariants;
    @NotNull
    private final List<SomaticVariant> consequencePassedVariants;

    VariantAnalysis(@NotNull final List<SomaticVariant> passedVariants,
            @NotNull final List<SomaticVariant> consensusPassedVariants,
            @NotNull final List<SomaticVariant> missenseVariants,
            @NotNull final List<SomaticVariant> consequencePassedVariants) {
        this.passedVariants = passedVariants;
        this.consensusPassedVariants = consensusPassedVariants;
        this.missenseVariants = missenseVariants;
        this.consequencePassedVariants = consequencePassedVariants;
    }

    @NotNull
    public List<SomaticVariant> passedVariants() {
        return passedVariants;
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
