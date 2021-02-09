package com.hartwig.hmftools.ckb.json;

import java.util.List;

import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.json.drug.Drug;
import com.hartwig.hmftools.ckb.json.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.json.gene.Gene;
import com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.reference.Reference;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;
import com.hartwig.hmftools.ckb.json.treatmentapproch.TreatmentApproach;
import com.hartwig.hmftools.ckb.json.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbJsonDatabase {

    @NotNull
    public abstract List<ClinicalTrial> clinicalTrial();

    @NotNull
    public abstract List<Drug> drug();

    @NotNull
    public abstract List<DrugClass> drugClass();

    @NotNull
    public abstract List<Gene> gene();

    @NotNull
    public abstract List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatus();

    @NotNull
    public abstract List<Indication> indication();

    @NotNull
    public abstract List<MolecularProfile> molecularProfile();

    @NotNull
    public abstract List<Reference> reference();

    @NotNull
    public abstract List<Therapy> therapy();

    @NotNull
    public abstract List<TreatmentApproach> treatmentApproach();

    @NotNull
    public abstract List<Variant> variant();
}
