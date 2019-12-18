package com.hartwig.hmftools.vicc.datamodel.civic;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicSource {

    @Nullable
    public abstract String name();

    @NotNull
    public abstract String status();

    @Nullable
    public abstract String openAccess();

    @Nullable
    public abstract String journal();

    @Nullable
    public abstract String fullJournalTitle();

    @NotNull
    public abstract String citation();

    @Nullable
    public abstract String pmcId();

    @NotNull
    public abstract String sourceUrl();

    @NotNull
    public abstract List<CivicClinicalTrial> clinicalTrials();

    @NotNull
    public abstract String pubmedId();

    @NotNull
    public abstract String isReview();

    @NotNull
    public abstract CivicPublicationDate publicationDate();

    @NotNull
    public abstract String id();
}
