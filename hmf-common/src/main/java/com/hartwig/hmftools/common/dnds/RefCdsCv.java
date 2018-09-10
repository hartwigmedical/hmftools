package com.hartwig.hmftools.common.dnds;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class RefCdsCv {

    public abstract String gene();

    public abstract int synonymousN();

    public abstract int missenseN();

    public abstract double missenseW();

    public abstract int nonsenseN();

    public abstract double nonsenseW();

    public abstract int spliceN();

    public abstract double spliceW();

    public abstract int indelN();

    public abstract double indelW();

    double missenseDrivers(int n) {
        return n * missenseProbability();
    }

    double missenseProbability() {
        return probability(missenseN(), missenseW());
    }

    private static double probability(int n, double w) {
        return n > 0 ? Math.max(0, (w - 1) / w) : 0;

    }

}
