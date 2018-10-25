package com.hartwig.hmftools.patientreporter.chord;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ChordAnalysis {

    public abstract double noneValue();

    public abstract double BRCA1Value();

    public abstract double BRCA2Value();

    public abstract double hrdValue();

    public abstract double predictedResponseValue();

}
