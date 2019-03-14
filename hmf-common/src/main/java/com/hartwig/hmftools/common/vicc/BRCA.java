package com.hartwig.hmftools.common.vicc;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BRCA {

    public abstract String variantFrequencyLOVD();

    public abstract String clinVarAccessionENIGMA();

    @Override
    public String toString() {
        return variantFrequencyLOVD() + ";" + clinVarAccessionENIGMA();
    }
}
