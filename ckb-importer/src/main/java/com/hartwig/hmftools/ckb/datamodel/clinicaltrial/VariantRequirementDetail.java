package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantRequirementDetail {

    public abstract int profileId();

    @NotNull
    public abstract String requirementType();

}