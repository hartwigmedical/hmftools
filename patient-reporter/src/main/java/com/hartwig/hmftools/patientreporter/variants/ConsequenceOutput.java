package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

class ConsequenceOutput {
    @NotNull
    private final List<SomaticVariant> consequentialVariants;
    @NotNull
    private final List<VariantReport> findings;

    ConsequenceOutput(@NotNull final List<SomaticVariant> consequentialVariants,
            @NotNull final List<VariantReport> findings) {
        this.consequentialVariants = consequentialVariants;
        this.findings = findings;
    }

    @NotNull
    List<SomaticVariant> consequentialVariants() {
        return consequentialVariants;
    }

    @NotNull
    List<VariantReport> findings() {
        return findings;
    }
}
