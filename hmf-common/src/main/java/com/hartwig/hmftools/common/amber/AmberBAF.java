package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.position.GenomePosition;

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
}
