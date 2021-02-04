package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.interpretation.clinicaltrialtree.ClinicalTrialInterpretation;
import com.hartwig.hmftools.ckb.interpretation.treatmenttree.TreatmentInterpretation;
import com.hartwig.hmftools.ckb.interpretation.varianttree.VariantInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntryInterpretation {

    @NotNull
    public abstract MolecularProfile molecularProfile();

    @NotNull
    public abstract List<VariantInterpretation> variantInterpretations();

    @NotNull
    public abstract List<TreatmentInterpretation> treatmentInterpretations();

    @NotNull
    public abstract List<ClinicalTrialInterpretation> clinicalTrialInterpretations();

}
