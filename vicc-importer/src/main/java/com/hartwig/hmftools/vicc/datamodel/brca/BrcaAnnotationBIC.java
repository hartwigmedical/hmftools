package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotationBIC {

    @NotNull
    public abstract String variantInBIC();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String mutationType();

    @NotNull
    public abstract String clinicalClassification();

    @NotNull
    public abstract String clinicalImportance();

    @NotNull
    public abstract String nomenclature();

    @NotNull
    public abstract String ethnicity();

    @NotNull
    public abstract String patientNationality();

    @NotNull
    public abstract String germlineOrSomatic();

    @NotNull
    public abstract String numberOfFamilyMemberCarryingMutation();

    @NotNull
    public abstract String literatureCitation();


}
