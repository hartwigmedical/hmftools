package com.hartwig.hmftools.lilac.util;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.ReferenceData.loadHlaTranscripts;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;

import java.util.HashMap;

import com.hartwig.hmftools.lilac.GeneCache;
import com.hartwig.hmftools.lilac.GeneSelector;

public class GeneCacheSetup
{
    public static void buildGeneCache()
    {
        if(GENE_CACHE != null)
            return;

        GENE_CACHE = new GeneCache(loadHlaTranscripts(V37, GeneSelector.MHC_CLASS_1));

        if(GENE_EXON_BOUNDARIES == null)
            GENE_EXON_BOUNDARIES = new HashMap<>();

        GENE_EXON_BOUNDARIES.put(HLA_A, GENE_CACHE.AminoAcidExonBoundaries.get(HLA_A));
        GENE_EXON_BOUNDARIES.put(HLA_B, GENE_CACHE.AminoAcidExonBoundaries.get(HLA_B));
        GENE_EXON_BOUNDARIES.put(HLA_C, GENE_CACHE.AminoAcidExonBoundaries.get(HLA_C));
    }
}
