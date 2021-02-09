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
import com.hartwig.hmftools.ckb.json.treatmentapproach.TreatmentApproach;
import com.hartwig.hmftools.ckb.json.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbJsonDatabase {

    @NotNull
    public abstract List<MolecularProfile> molecularProfiles();

    @NotNull
    public abstract List<Variant> variants();

    @NotNull
    public abstract List<Gene> genes();

    @NotNull
    public abstract List<Indication> indications();

    @NotNull
    public abstract List<TreatmentApproach> treatmentApproaches();

    @NotNull
    public abstract List<Therapy> therapies();

    @NotNull
    public abstract List<Drug> drugs();

    @NotNull
    public abstract List<DrugClass> drugClasses();

    @NotNull
    public abstract List<ClinicalTrial> clinicalTrials();

    @NotNull
    public abstract List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatuses();

    @NotNull
    public abstract List<Reference> references();
}
