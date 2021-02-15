package com.hartwig.hmftools.ckb.json;

import java.util.List;

import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrial;
import com.hartwig.hmftools.ckb.json.drug.JsonDrug;
import com.hartwig.hmftools.ckb.json.drugclass.JsonDrugClass;
import com.hartwig.hmftools.ckb.json.gene.JsonGene;
import com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus.JsonGlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.json.indication.JsonIndication;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;
import com.hartwig.hmftools.ckb.json.reference.JsonReference;
import com.hartwig.hmftools.ckb.json.therapy.JsonTherapy;
import com.hartwig.hmftools.ckb.json.treatmentapproach.JsonTreatmentApproach;
import com.hartwig.hmftools.ckb.json.variant.JsonVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbJsonDatabase {

    @NotNull
    public abstract List<JsonMolecularProfile> molecularProfiles();

    @NotNull
    public abstract List<JsonVariant> variants();

    @NotNull
    public abstract List<JsonGene> genes();

    @NotNull
    public abstract List<JsonIndication> indications();

    @NotNull
    public abstract List<JsonTreatmentApproach> treatmentApproaches();

    @NotNull
    public abstract List<JsonTherapy> therapies();

    @NotNull
    public abstract List<JsonDrug> drugs();

    @NotNull
    public abstract List<JsonDrugClass> drugClasses();

    @NotNull
    public abstract List<JsonClinicalTrial> clinicalTrials();

    @NotNull
    public abstract List<JsonGlobalTherapyApprovalStatus> globalTherapyApprovalStatuses();

    @NotNull
    public abstract List<JsonReference> references();
}
