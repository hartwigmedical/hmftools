package com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileInterpretation {

    @Nullable
    public abstract List<Variant> variants();
}
