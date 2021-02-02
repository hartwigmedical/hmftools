package com.hartwig.hmftools.ckb.interpretation;

import com.hartwig.hmftools.ckb.datamodel.molecularprofile.MolecularProfile;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntryInterpretation {

    @NotNull
    public abstract MolecularProfile molecularProfile();

    @NotNull
    public abstract VariantInterpretation variantInterpretation();

}
