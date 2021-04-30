package com.hartwig.hmftools.protect.viralbreakend;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VirusBreakend {

    @NotNull
    public abstract String taxidGenus();

    @NotNull
    public abstract String nameGenus();

    @NotNull
    public abstract String readsGenusTree();

    @NotNull
    public abstract String taxidSpecies();

    @NotNull
    public abstract String nameSpecies();

    @NotNull
    public abstract String readsSpeciesTree();

    @NotNull
    public abstract String taxidAssigned();

    @NotNull
    public abstract String nameAssigned();

    @NotNull
    public abstract String readsAssignedTree();

    @NotNull
    public abstract String readsAssignedDirect();

    @NotNull
    public abstract String Reference();

    @NotNull
    public abstract String referenceTaxid();

    @NotNull
    public abstract String referenceKmerCount();

    @NotNull
    public abstract String alternateKmerCountRname();

    @NotNull
    public abstract String startpos();

    @NotNull
    public abstract String endpos();

    @NotNull
    public abstract String numreads();

    @NotNull
    public abstract String covbases();

    @NotNull
    public abstract String coverage();

    @NotNull
    public abstract String meandepth();

    @NotNull
    public abstract String meanbaseq();

    @NotNull
    public abstract String meanmapq();

    @NotNull
    public abstract String integrations();

    @NotNull
    public abstract String QCStatus();
}