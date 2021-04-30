package com.hartwig.hmftools.common.virusbreakend;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VirusBreakend {

    public abstract int taxidGenus();

    @NotNull
    public abstract String nameGenus();

    public abstract int readsGenusTree();

    public abstract int taxidSpecies();

    @NotNull
    public abstract String nameSpecies();

    public abstract int readsSpeciesTree();

    public abstract int taxidAssigned();

    @NotNull
    public abstract String nameAssigned();

    public abstract int readsAssignedTree();

    public abstract int readsAssignedDirect();

    @NotNull
    public abstract String reference();

    public abstract int referenceTaxid();

    public abstract int referenceKmerCount();

    @NotNull
    public abstract String alternateKmerCountRname();

    public abstract int startpos();

    public abstract int endpos();

    public abstract int numreads();

    public abstract int covbases();
    
    public abstract int coverage();

    public abstract int meandepth();

    public abstract int meanbaseq();

    public abstract int meanmapq();

    public abstract int integrations();

    public abstract int QCStatus();
}