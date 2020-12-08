package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneCheckerTestFactory;

import org.junit.Test;

public class FusionExtractorTest {

    private static final GeneChecker HG19_GENE_CHECKER = GeneCheckerTestFactory.buildForHG19();

    @Test
    public void canExtractFusionPairsGenes() {
        FusionExtractor fusionExtractor = new FusionExtractor(HG19_GENE_CHECKER);
        KnownFusionPair fusion = fusionExtractor.extract(Knowledgebase.VICC_CGI, "PDGFRA", EventType.FUSION_PAIR, "BCR-PDGFRA Fusion");

        assertNotNull(fusion);
        assertEquals("BCR", fusion.geneUp());
        assertEquals("PDGFRA", fusion.geneDown());
    }

    @Test
    public void ignoresFusionsOnUnknownGenes() {
        FusionExtractor fusionExtractor = new FusionExtractor(HG19_GENE_CHECKER);
        KnownFusionPair fusion = fusionExtractor.extract(Knowledgebase.VICC_CGI, "IG", EventType.FUSION_PAIR, "IG-BCL2");

        assertNull(fusion);
    }

    @Test
    public void canExtractFusionPairsWithExonsUpDown() {
        FusionExtractor fusionExtractor = new FusionExtractor(HG19_GENE_CHECKER);
        KnownFusionPair fusion = fusionExtractor.extract(Knowledgebase.VICC_CGI, "EGFR", EventType.FUSION_PAIR, "EGFRvII");

        assertNotNull(fusion);
        assertEquals("EGFR", fusion.geneUp());
        assertEquals(13, (int) fusion.minExonUp());
        assertEquals(13, (int) fusion.maxExonUp());
        assertEquals("EGFR", fusion.geneDown());
        assertEquals(16, (int) fusion.minExonDown());
        assertEquals(16, (int) fusion.maxExonDown());
    }

    @Test
    public void canExtractFusionPairsWithOddNames() {
        FusionExtractor fusionExtractor = new FusionExtractor(new GeneChecker(Sets.newHashSet("IGH", "NKX2-1")));
        KnownFusionPair fusion = fusionExtractor.extract(Knowledgebase.VICC_CGI, "NKX2-1", EventType.FUSION_PAIR, "IGH-NKX2-1 Fusion");

        assertNotNull(fusion);
        assertEquals("IGH", fusion.geneUp());
        assertEquals("NKX2-1", fusion.geneDown());
    }

    @Test
    public void canExtractFusionPairsWithRangeExons() {
        FusionExtractor fusionExtractor = new FusionExtractor(HG19_GENE_CHECKER);
        KnownFusionPair fusion =
                fusionExtractor.extract(Knowledgebase.VICC_CGI, "MET", EventType.FUSION_PAIR_AND_EXON, "EXON 14 SKIPPING MUTATION");

        assertNotNull(fusion);
        assertEquals("MET", fusion.geneUp());
        assertEquals(13, (int) fusion.minExonUp());
        assertEquals(13, (int) fusion.maxExonUp());
        assertEquals("MET", fusion.geneDown());
        assertEquals(15, (int) fusion.minExonDown());
        assertEquals(15, (int) fusion.maxExonDown());
    }
}