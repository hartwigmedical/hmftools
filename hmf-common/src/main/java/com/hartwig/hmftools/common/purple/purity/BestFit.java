package com.hartwig.hmftools.common.purple.purity;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface BestFit {

    @NotNull
    FittedPurity fit();

    @NotNull
    FittedPurityScore score();

    @NotNull
    FittedPurityStatus status();

    @NotNull
    List<FittedPurity> bestFitPerPurity();
}
