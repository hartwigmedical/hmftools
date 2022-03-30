package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AmberBAF extends GenomePosition {

    double tumorBAF();

    default double tumorModifiedBAF() {
        return 0.5 + Math.abs(tumorBAF() - 0.5);
    }

    int tumorDepth();

    double normalBAF();

    default double normalModifiedBAF() {
        return 0.5 + Math.abs(normalBAF() - 0.5);
    }

    int normalDepth();

    @NotNull
    static AmberBAF create(@NotNull final TumorBAF tumor) {
        int tumorAltCount = tumor.tumorAltSupport();
        double tumorBaf = tumorAltCount / (double) (tumorAltCount + tumor.tumorRefSupport());
        int normalAltCount = tumor.normalAltSupport();
        double normalBaf = normalAltCount / (double) (normalAltCount + tumor.normalRefSupport());
        return ImmutableAmberBAF.builder()
                .from(tumor)
                .normalDepth(tumor.normalReadDepth())
                .tumorDepth(tumor.tumorReadDepth())
                .normalBAF(normalBaf)
                .tumorBAF(tumorBaf)
                .build();
    }

    @NotNull
    static AmberBAF create(@NotNull final BaseDepth baseDepth) {
        int normalAltCount = baseDepth.altSupport();
        double normalBaf = normalAltCount / (double) (normalAltCount + baseDepth.refSupport());
        return ImmutableAmberBAF.builder()
                .from(baseDepth)
                .normalDepth(baseDepth.readDepth())
                .tumorDepth(-1)
                .normalBAF(normalBaf)
                .tumorBAF(-1)
                .build();
    }
}
