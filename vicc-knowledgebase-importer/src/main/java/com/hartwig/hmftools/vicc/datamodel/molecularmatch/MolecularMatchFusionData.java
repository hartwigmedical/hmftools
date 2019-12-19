package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchFusionData {

    @Nullable
    public abstract String source();

    @Nullable
    public abstract String synonym();

    @NotNull
    public abstract List<String> aChromosomes();

    @NotNull
    public abstract List<String> aBands();

    @NotNull
    public abstract List<String> aGenes();

    @NotNull
    public abstract List<String> aCoords();

    @NotNull
    public abstract List<String> aTranscripts();

    @NotNull
    public abstract List<String> aOrientations();

    @NotNull
    public abstract List<MolecularMatchFusionGenomicRegion> aGenomicRegions();

    @NotNull
    public abstract List<String> bChromosomes();

    @NotNull
    public abstract List<String> bBands();

    @NotNull
    public abstract List<String> bGenes();

    @NotNull
    public abstract List<String> bCoords();

    @NotNull
    public abstract List<String> bTranscripts();

    @NotNull
    public abstract List<String> bOrientations();

    @NotNull
    public abstract List<MolecularMatchFusionGenomicRegion> bGenomicRegions();

    @NotNull
    public abstract List<String> inserts();

    @Nullable
    public abstract String paper();
}
