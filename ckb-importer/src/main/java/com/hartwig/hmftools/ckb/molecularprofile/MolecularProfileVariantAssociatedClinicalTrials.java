package com.hartwig.hmftools.ckb.molecularprofile;

import java.util.List;

import com.hartwig.hmftools.ckb.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileVariantAssociatedClinicalTrials {

    @NotNull
    public abstract String nctId();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract String phase();

    @NotNull
    public abstract String recruitment();

    @NotNull
    public abstract List<TherapyInfo> therapy();

}
