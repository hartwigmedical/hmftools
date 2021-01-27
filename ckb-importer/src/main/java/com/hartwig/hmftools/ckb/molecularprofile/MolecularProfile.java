package com.hartwig.hmftools.ckb.molecularprofile;

import java.util.List;

import com.hartwig.hmftools.ckb.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfile {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract List<VariantInfo> geneVariant();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproach();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract String updateDate();

    @NotNull
    public abstract MolecularProfileExtendedEvidence complexMolecularProfileEvidence();

    @NotNull
    public abstract MolecularProfileExtendedEvidence treatmentApproachEvidence();

    @NotNull
    public abstract List<ClinicalTrialInfo> variantAssociatedClinicalTrial();

    @NotNull
    public abstract MolecularProfileExtendedEvidence variantLevelEvidence();

    @NotNull
    public abstract MolecularProfileExtendedEvidence extendedEvidence();

}
