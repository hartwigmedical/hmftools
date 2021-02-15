package com.hartwig.hmftools.ckb.datamodel.common.molecularprofile;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileInterpretation {

    @Nullable
    public abstract List<Variant> variants();
}
