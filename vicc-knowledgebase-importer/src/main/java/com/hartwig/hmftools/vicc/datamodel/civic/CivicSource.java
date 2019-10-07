package com.hartwig.hmftools.vicc.datamodel.civic;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicSource {

    @NotNull
    public abstract String status();

    @Nullable
    public abstract String openAccess();

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String journal();

    @NotNull
    public abstract String citation();

    @Nullable
    public abstract String pmc_Id();

    @Nullable
    public abstract String fullJournalTitle();

    @NotNull
    public abstract String sourceUrl();

    @Nullable
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
