package com.hartwig.hmftools.ckb.json.clinicaltrial;

import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrialVariantRequirementDetail {

    @NotNull
    public abstract MolecularProfileInfo molecularProfile();

    @NotNull
    public abstract String requirementType();
}
