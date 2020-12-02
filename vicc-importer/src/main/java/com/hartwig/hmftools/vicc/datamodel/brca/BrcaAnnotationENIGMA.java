package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotationENIGMA {

    @NotNull
    public abstract String variantInENIGMA();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String alleleOrigin();

    @NotNull
    public abstract String clinVarAccession();

    @NotNull
    public abstract String assertionMethod();

    @NotNull
    public abstract String assertionMethodCitation();

    @NotNull
    public abstract String collectionMethod();

    @NotNull
    public abstract String conditionCategory();

    @NotNull
    public abstract String conditionIdValue();

    @NotNull
    public abstract String conditionIdType();

    @NotNull
    public abstract String clinicalSignificance();

    @NotNull
    public abstract String clinicalSignificanceCitations();

    @NotNull
    public abstract String commentOnClinicalSignificance();

    @NotNull
    public abstract String dateLastEvaluated();

    @NotNull
    public abstract String url();


}
