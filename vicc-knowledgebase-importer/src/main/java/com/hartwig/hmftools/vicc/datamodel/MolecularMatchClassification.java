package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchClassification {

    @Nullable
    public abstract List<String> end();

    @NotNull
    public abstract String classification();

    @Nullable
    public abstract String classificationOverride();

    @Nullable
    public abstract List<String> start();

    @Nullable
    public abstract List<String> chr();

    @Nullable
    public abstract String geneSymbol();

    @Nullable
    public abstract List<String> pathology();

    @Nullable
    public abstract List<String> ref();

    @Nullable
    public abstract String description();

    @Nullable
    public abstract String priority();

    @Nullable
    public abstract List<String> NucleotideChange();

    @Nullable
    public abstract String expandGeneSearch();

    @Nullable
    public abstract List<MolecularMatchParents> parents();

    @Nullable
    public abstract String drugsExperimentalCount();

    @Nullable
    public abstract List<String> exon();

    @Nullable
    public abstract String drugsApprovedOffLabelCount();

    @Nullable
    public abstract List<String> exonicFunc();

    @Nullable
    public abstract List<String> popFreqMax();

    @Nullable
    public abstract String copyNumberType();

    @Nullable
    public abstract String publicationCount();

    @Nullable
    public abstract String transcript();

    @Nullable
    public abstract List<String> dbSNP();

    @Nullable
    public abstract List<String> alt();

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String rootTerm();

    @Nullable
    public abstract List<String> sources();

    @Nullable
    public abstract String drugsApprovedOnLabelCount();

    @Nullable
    public abstract String trialCount();

    @Nullable
    public abstract String alias();

    @Nullable
    public abstract List<String> COSMIC_ID();

    @Nullable
    public abstract List<String> transcripts();
}
