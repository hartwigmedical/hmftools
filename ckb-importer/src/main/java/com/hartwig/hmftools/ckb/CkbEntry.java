package com.hartwig.hmftools.ckb;

import java.util.List;

import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.drug.Drug;
import com.hartwig.hmftools.ckb.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.gene.Gene;
import com.hartwig.hmftools.ckb.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.indication.Indication;
import com.hartwig.hmftools.ckb.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.reference.Reference;
import com.hartwig.hmftools.ckb.therapy.Therapy;
import com.hartwig.hmftools.ckb.treatmentApproach.TreatmentApproach;
import com.hartwig.hmftools.ckb.variant.Variant;

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
