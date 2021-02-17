package com.hartwig.hmftools.ckb.json.molecularprofile;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JsonMolecularProfile implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract List<VariantInfo> geneVariants();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproaches();

    @Nullable
    public abstract LocalDate createDate();

    @Nullable
    public abstract LocalDate updateDate();

    @NotNull
    public abstract JsonMolecularProfileExtendedEvidence complexMolecularProfileEvidence();

    @NotNull
    public abstract JsonMolecularProfileExtendedEvidence treatmentApproachEvidence();

    @NotNull
    public abstract List<ClinicalTrialInfo> variantAssociatedClinicalTrials();

    @NotNull
    public abstract JsonMolecularProfileExtendedEvidence variantLevelEvidence();

    @NotNull
    public abstract JsonMolecularProfileExtendedEvidence extendedEvidence();
}
