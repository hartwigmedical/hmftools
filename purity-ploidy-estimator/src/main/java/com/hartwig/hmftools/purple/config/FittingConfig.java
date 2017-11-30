package com.hartwig.hmftools.purple.config;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FittingConfig {

    double minPurity();

    double maxPurity();

    double purityIncrement();

    double minNormFactor();

    double maxNormFactor();

    double normFactorIncrement();

    default int maxPloidy() {
        return 20;
    }
}
