package com.hartwig.hmftools.ckb.datamodelinterpretation;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.datamodelinterpretation.gene.Gene;
import com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.Indication;
import com.hartwig.hmftools.ckb.datamodelinterpretation.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.datamodelinterpretation.reference.Reference;
import com.hartwig.hmftools.ckb.datamodelinterpretation.treatmentApproch.TreatmentApproch;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntryInterpretation {

    public abstract int id();

    @NotNull
    public abstract Gene gene();

    @NotNull
    public abstract List<Variant> variants(); //complex is possible

    @NotNull
    public abstract Therapy therapy();

    @NotNull
    public abstract List<Drug> drugs();

    @NotNull
    public abstract DrugClass drugClass();

    @NotNull
    public abstract List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatus();

    @NotNull
    public abstract TreatmentApproch treatmentApproch();

    @NotNull
    public abstract ClinicalTrial clinicalTrial();

    @NotNull
    public abstract Indication indication();

    @NotNull
    public abstract MolecularProfile molecularProfile();

    @NotNull
    public abstract List<Reference> reference();

}
