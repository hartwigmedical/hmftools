package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchClassification {

    @NotNull
    public abstract List<String> end();

    @NotNull
    public abstract String classification();

    @NotNull
    public abstract String classificationOverride();

    @NotNull
    public abstract List<String> start();

    @NotNull
    public abstract List<String> chr();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> pathology();

    @NotNull
    public abstract List<String> ref();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String priority();

    @NotNull
    public abstract List<String> NucleotideChange();

    @NotNull
    public abstract List<String> parents();

    @NotNull
    public abstract String drugsExperimentalCount();

    @NotNull
    public abstract List<String> exon();

    @NotNull
    public abstract String drugsApprovedOffLabelCount();

    @NotNull
    public abstract List<String> exonicFunc();

    @NotNull
    public abstract List<String> popFreqMax();

    @NotNull
    public abstract String copyNumberType();

    @NotNull
    public abstract String publicationCount();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract List<String> dbSNP();

    @NotNull
    public abstract List<String> alt();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String rootTerm();

    @NotNull
    public abstract List<String> sources();

    @NotNull
    public abstract String drugsApprovedOnLabelCount();

    @NotNull
    public abstract String trialCount();

    @NotNull
    public abstract String alias();

    @NotNull
    public abstract List<String> COSMIC_ID();

    @NotNull
    public abstract List<String> transcripts();
}
