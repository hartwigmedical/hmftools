package com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.Variant;
import com.hartwig.hmftools.ckb.interpretation.VariantInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrialVariantRequirementDetail {

    public abstract int id();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract String requirementType();

    @Nullable
    public abstract List<VariantInterpretation> variantInterpretation();
}
