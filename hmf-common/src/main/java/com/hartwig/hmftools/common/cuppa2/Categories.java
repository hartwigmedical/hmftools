package com.hartwig.hmftools.common.cuppa2;

public class Categories {

    public enum ClfName {
        combined,

        dna_combined,
        gen_pos,
        snv96,
        event,

        rna_combined,
        gene_exp,
        alt_sj,

        none
        ;
    }

    public enum ClfGroup {
        combined,
        dna,
        rna,

        none
        ;
    }

    public enum DataType {
        prob,
        sig_quantile,
        feat_contrib,
        cv_performance,

        none
        ;
    }
}
