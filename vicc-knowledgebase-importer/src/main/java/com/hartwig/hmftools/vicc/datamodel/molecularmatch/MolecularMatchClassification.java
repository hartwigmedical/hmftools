package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchClassification {

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String geneSymbol();

    @Nullable
    public abstract String expandGeneSearch();

    @Nullable
    public abstract String transcript();

    @NotNull
    public abstract List<String> transcripts();

    @NotNull
    public abstract List<String> chromosomes();

    @NotNull
    public abstract List<String> starts();

    @NotNull
    public abstract List<String> ends();

    @NotNull
    public abstract List<String> refs();

    @NotNull
    public abstract List<String> alts();

    @NotNull
    public abstract List<String> nucleotideChanges();

    @NotNull
    public abstract List<String> exons();

    @NotNull
    public abstract List<String> exonicFuncs();

    @NotNull
    public abstract String classification();

    @Nullable
    public abstract String classificationOverride();

    @NotNull
    public abstract List<String> pathology();

    @Nullable
    public abstract String copyNumberType();

    @Nullable
    public abstract String drugsApprovedOnLabelCount();

    @Nullable
    public abstract String drugsApprovedOffLabelCount();

    @Nullable
    public abstract String drugsExperimentalCount();

    @Nullable
    public abstract String trialCount();

    @Nullable
    public abstract String publicationCount();

    @NotNull
    public abstract List<String> sources();

    @NotNull
    public abstract List<String> dbSNPs();

    @NotNull
    public abstract List<String> cosmicIds();

    @NotNull
    public abstract List<String> popFreqMaxes();

    @NotNull
    public abstract List<MolecularMatchParent> parents();

    @Nullable
    public abstract String rootTerm();

    @Nullable
    public abstract String alias();

    @Nullable
    public abstract String priority();

    @Nullable
    public abstract String description();
}
