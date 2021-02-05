package com.hartwig.hmftools.ckb.datamodel;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.gene.Gene;
import com.hartwig.hmftools.ckb.datamodel.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.datamodel.indication.Indication;
import com.hartwig.hmftools.ckb.datamodel.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.TreatmentApproach;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntry {

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
