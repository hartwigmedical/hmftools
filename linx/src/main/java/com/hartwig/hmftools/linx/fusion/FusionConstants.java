package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_LINC_RNA;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROCESSED_TRANS;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_RETAINED_INTRON;

import java.util.List;

import com.google.common.collect.Lists;

public class FusionConstants
{
    // the maximum distance upstream of a gene for a breakend to be consider a fusion candidate
    public static final int DEFAULT_PRE_GENE_PROMOTOR_DISTANCE = 100000;

    public static int PRE_GENE_PROMOTOR_DISTANCE = DEFAULT_PRE_GENE_PROMOTOR_DISTANCE; // can be overridden in config

    // maximum distance from breakend to transcript start for known and other fusions
    public static final int MAX_UPSTREAM_DISTANCE_KNOWN = 100000;
    public static final int MAX_UPSTREAM_DISTANCE_OTHER = 10000;

    public static final int SHORT_UNPHASED_DISTANCE_KNOWN = 1000000;

    public static final List<String> REQUIRED_BIOTYPES = Lists.newArrayList(
            BIOTYPE_PROCESSED_TRANS, BIOTYPE_PROTEIN_CODING, BIOTYPE_NONSENSE_MED_DECAY, BIOTYPE_RETAINED_INTRON, BIOTYPE_LINC_RNA);

    public static final int FUSION_MAX_CHAIN_LENGTH = 150000;
    public static final int FUSION_MAX_CHAIN_LINKS = 4;

}
