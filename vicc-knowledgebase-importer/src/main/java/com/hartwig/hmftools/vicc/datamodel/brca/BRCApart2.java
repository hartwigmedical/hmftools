package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BRCApart2 {

    @NotNull
    public abstract String Allele_count_FIN_ExAC();

    @NotNull
    public abstract String Condition_category_ENIGMA();

    @NotNull
    public abstract String Allele_Frequency_ESP();

    @NotNull
    public abstract String Homozygous_count_OTH_ExAC();

    @NotNull
    public abstract String Genetic_origin_LOVD();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String Homozygous_count_AMR_ExAC();

    @NotNull
    public abstract String Clinical_Significance_ClinVar();

    @NotNull
    public abstract String AA_Allele_Frequency_ESP();

    @NotNull
    public abstract String Protein_Change();

    @NotNull
    public abstract String Variant_in_exLOVD();

    @NotNull
    public abstract String EA_Allele_Frequency_ESP();

    @NotNull
    public abstract String HGVS_RNA();

    @NotNull
    public abstract String Clinical_significance_citations_ENIGMA();

    @NotNull
    public abstract String Variant_effect_LOVD();

    @NotNull
    public abstract String Polyphen_Prediction();

    @NotNull
    public abstract String Data_Release_id();

    @NotNull
    public abstract String Hg37_Start();

    @NotNull
    public abstract String Hg36_Start();

    @NotNull
    public abstract String Sift_Score();

    @NotNull
    public abstract String Genomic_Coordinate_hg38();

    @NotNull
    public abstract String Alt();

    @NotNull
    public abstract String Literature_citation_BIC();

    @NotNull
    public abstract String Variant_haplotype_LOVD();

    @NotNull
    public abstract String Allele_frequency_NFE_ExAC();

    @NotNull
    public abstract String Hg38_Start();

    @NotNull
    public abstract String Pos();

    @NotNull
    public abstract String Date_last_evaluated_ENIGMA();

    @NotNull
    public abstract String Allele_number_SAS_ExAC();

    @NotNull
    public abstract String Allele_number_AMR_ExAC();

    @NotNull
    public abstract String DBID_LOVD();
}
