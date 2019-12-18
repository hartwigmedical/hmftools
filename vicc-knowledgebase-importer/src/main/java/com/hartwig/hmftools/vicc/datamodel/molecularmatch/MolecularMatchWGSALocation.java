package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchWGSALocation {

    @NotNull
    public abstract List<String> genes();

    @NotNull
    public abstract String chr();

    @NotNull
    public abstract String start();

    @NotNull
    public abstract String end();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String chrStartRefAlt();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String nucleotideChange();

    @Nullable
    public abstract String aa();

    @NotNull
    public abstract List<String> fullAAs();

    @Nullable
    public abstract String exonicFunc();

    @NotNull
    public abstract String popFreqMax();

    @NotNull
    public abstract List<String> clinVarDiseases();

    @NotNull
    public abstract List<String> clinVarSigs();

    @NotNull
    public abstract List<String> clinVarStates();

    @NotNull
    public abstract List<String> clinVarDbIds();

    @Nullable
    public abstract String exacAFR();

    @Nullable
    public abstract String exacAMR();

    @Nullable
    public abstract String exacEAS();

    @Nullable
    public abstract String exacFIN();

    @Nullable
    public abstract String exacNFE();

    @Nullable
    public abstract String exacSAS();

    @Nullable
    public abstract String exacFreq();

    @Nullable
    public abstract String g1000AFR();

    @Nullable
    public abstract String g1000AMR();

    @Nullable
    public abstract String g1000EAS();

    @Nullable
    public abstract String g1000EUR();

    @Nullable
    public abstract String g1000SAS();

    @Nullable
    public abstract String g1000ALL();

    @NotNull
    public abstract String fathmm();

    @NotNull
    public abstract String fathmmPred();

    @Nullable
    public abstract String esp6500siAA();

    @Nullable
    public abstract String esp6500siEA();

    @Nullable
    public abstract String dbSNP();

    @Nullable
    public abstract String cosmicId();

    @NotNull
    public abstract String phyloP46wayPlacental();

    @NotNull
    public abstract String phyloP100wayVertebrate();

    @NotNull
    public abstract String siPhy29wayLogOdds();

    @Nullable
    public abstract String gwasSNP();

    @Nullable
    public abstract String gwasDIS();

    @Nullable
    public abstract String gwasPubmed();

    @NotNull
    public abstract String gerpRS();

    @NotNull
    public abstract String func();

    @Nullable
    public abstract String wgRna();

    @Nullable
    public abstract String targetScanS();

    @NotNull
    public abstract String key();

}
