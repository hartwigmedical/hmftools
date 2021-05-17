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

    public abstract int alternateKmerCount();

    @NotNull
    public abstract String RName();

    public abstract int startPos();

    public abstract int endPos();

    public abstract int numReads();

    public abstract int covBases();
    
    public abstract double coverage();

    public abstract double meanDepth();

    public abstract double meanBaseQ();

    public abstract double meanMapQ();

    public abstract int integrations();

    @NotNull
    public abstract VirusBreakendQCStatus qcStatus();
}