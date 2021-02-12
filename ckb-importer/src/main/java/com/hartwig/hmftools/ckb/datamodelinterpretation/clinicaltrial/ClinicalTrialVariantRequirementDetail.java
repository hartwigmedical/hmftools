package com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial;

import com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation.MolecularProfileInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrialVariantRequirementDetail {

    public abstract int id();  // TODO: could this removed, because is this present in MolecularProfileInterpretation object?

    @NotNull
    public abstract String profileName(); // TODO: could this removed, because is this present in MolecularProfileInterpretation object?

    @NotNull
    public abstract String requirementType();

    @Nullable
    public abstract MolecularProfileInterpretation molecularProfileInterpretation();
}