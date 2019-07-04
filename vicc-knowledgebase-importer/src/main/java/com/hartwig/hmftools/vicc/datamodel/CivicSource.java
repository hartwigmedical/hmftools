package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicSource {

    @NotNull
    public abstract String status();

    @NotNull
    public abstract String openAccess();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String journal();

    @NotNull
    public abstract String citation();

    @NotNull
    public abstract String pmc_I();

    @NotNull
    public abstract String fullJournalTitle();

    @NotNull
    public abstract String sourceUrl();

    @NotNull
    public abstract List<String> clinicalTrials();

    @NotNull
    public abstract String pubmedId();

    @NotNull
    public abstract String isReview();

    @NotNull
    public abstract CivicPublicationDate publicationDate();

    @NotNull
    public abstract String id();
}
