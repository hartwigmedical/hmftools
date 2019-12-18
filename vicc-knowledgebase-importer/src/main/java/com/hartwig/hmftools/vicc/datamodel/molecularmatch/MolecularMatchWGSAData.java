package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchWGSAData {

    @Nullable
    public abstract String ExonicFunc();

    @Nullable
    public abstract String dbSNP();

    @Nullable
    public abstract List<String> ClinVar_DIS();

    @Nullable
    public abstract List<String> ClinVar_SIG();

    @Nullable
    public abstract List<String> ClinVar_STATUS();

    @Nullable
    public abstract List<String> ClinVar_DBID();

    @Nullable
    public abstract String ExAC_NFE();

    @Nullable
    public abstract String ExAC_FIN();

    @Nullable
    public abstract String G1000_ALL();

    @Nullable
    public abstract String G1000_SAS();

    @Nullable
    public abstract String G1000_EAS();

    @Nullable
    public abstract String G1000_AFR();

    @Nullable
    public abstract String ExAC_SAS();

    @Nullable
    public abstract String ExAC_EAS();

    @Nullable
    public abstract String ExAC_AMR();

    @Nullable
    public abstract String ExAC_AFR();

    @Nullable
    public abstract String ExAC_Freq();

    @NotNull
    public abstract String End();

    @NotNull
    public abstract String Start();

    @NotNull
    public abstract String SiPhy_29way_logOdds();

    @NotNull
    public abstract List<String> FullAA();

    @NotNull
    public abstract String Ref();

    @NotNull
    public abstract String GERP_RS();

    @NotNull
    public abstract String FATHMM();

    @NotNull
    public abstract String NucleotideChange();

    @NotNull
    public abstract String phyloP100way_vertebrate();

    @NotNull
    public abstract String Func();

    @Nullable
    public abstract String GWAS_PUBMED();

    @NotNull
    public abstract String Transcript();

    @Nullable
    public abstract String ESP6500si_AA();

    @Nullable
    public abstract String ESP6500si_EA();

    @Nullable
    public abstract String G1000_EUR();

    @Nullable
    public abstract String G1000_AMR();

    @NotNull
    public abstract String Chr_Start_Ref_Alt();

    @Nullable
    public abstract String AA();

    @NotNull
    public abstract String PopFreqMax();

    @NotNull
    public abstract String FATHMM_Pred();

    @Nullable
    public abstract String wgRna();

    @NotNull
    public abstract List<String> Gene();

    @NotNull
    public abstract String phyloP46way_placental();

    @NotNull
    public abstract String key();

    @Nullable
    public abstract String targetScanS();

    @NotNull
    public abstract String Chr();

    @Nullable
    public abstract String COSMIC_ID();

    @NotNull
    public abstract String alt();

    @Nullable
    public abstract String GWAS_DIS();

    @Nullable
    public abstract String GWAS_SNP();

}
