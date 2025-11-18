package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_LINC_RNA;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROCESSED_TRANS;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_RETAINED_INTRON;
import static com.hartwig.hmftools.common.fusion.FusionCommon.DEFAULT_PRE_GENE_PROMOTOR_DISTANCE;

import java.util.List;

import com.google.common.collect.Lists;

public class FusionConstants
{
    // the maximum distance upstream of a gene for a breakend to be considered a fusion candidate

    public static int PRE_GENE_PROMOTOR_DISTANCE = DEFAULT_PRE_GENE_PROMOTOR_DISTANCE; // can be overridden in config

    // maximum distance from breakend to transcript start for known and other fusions
    public static final int MAX_UPSTREAM_DISTANCE_KNOWN = 100000;
    public static final int MAX_UPSTREAM_DISTANCE_IG_KNOWN = 250000;
    public static final int MAX_UPSTREAM_DISTANCE_OTHER = 10000;

    public static final int SHORT_UNPHASED_DISTANCE_KNOWN = 1000000;

    public static final List<String> REQUIRED_BIOTYPES = Lists.newArrayList(
            BIOTYPE_PROCESSED_TRANS, BIOTYPE_PROTEIN_CODING, BIOTYPE_NONSENSE_MED_DECAY, BIOTYPE_RETAINED_INTRON, BIOTYPE_LINC_RNA);

    public static final int FUSION_MAX_CHAIN_LENGTH = 150000;
    public static final int FUSION_MAX_CHAIN_LINKS = 4;

    public static final int ENHANCER_PROMISCUOUS_MIN_DISTANCE = 100000;

    public static final List<String> PROTEINS_REQUIRED_KEPT = Lists.newArrayList();

    public static final List<String> PROTEINS_REQUIRED_LOST = Lists.newArrayList();

    static
    {
        // first entries are for V37, second are for V38 in Ensembl
        PROTEINS_REQUIRED_LOST.add("Raf-like Ras-binding");
        PROTEINS_REQUIRED_LOST.add("RBD");

        PROTEINS_REQUIRED_KEPT.add("Ets domain");
        PROTEINS_REQUIRED_KEPT.add("ETS_DOMAIN_3");

        PROTEINS_REQUIRED_KEPT.add("Protein kinase domain");
        PROTEINS_REQUIRED_KEPT.add("PROTEIN_KINASE_DOM");

        PROTEINS_REQUIRED_KEPT.add("Epidermal growth factor-like domain");
        PROTEINS_REQUIRED_KEPT.add("EGF_3");

        PROTEINS_REQUIRED_KEPT.add("Ankyrin repeat-containing domain");
        PROTEINS_REQUIRED_KEPT.add("ANK_REPEAT");
        PROTEINS_REQUIRED_KEPT.add("ANK_REP_REGION");

        PROTEINS_REQUIRED_KEPT.add("Basic-leucine zipper domain");
        PROTEINS_REQUIRED_KEPT.add("BZIP");

        PROTEINS_REQUIRED_KEPT.add("High mobility group box domain");
        PROTEINS_REQUIRED_KEPT.add("HMG_BOX_2");

        PROTEINS_REQUIRED_KEPT.add("Bromodomain");
        PROTEINS_REQUIRED_KEPT.add("BROMODOMAIN_2");
    }
}
