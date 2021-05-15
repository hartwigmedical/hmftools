package com.hartwig.hmftools.serve.extraction.fusion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.extraction.util.GeneCheckerTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FusionExtractorTest {

    private static final GeneChecker V37_GENE_CHECKER = GeneCheckerTestFactory.buildForV37();

    @Test
    public void canExtractSimpleFusionPair() {
        FusionExtractor fusionExtractor = testFusionExtractor();
        KnownFusionPair fusion = fusionExtractor.extract("PDGFRA", EventType.FUSION_PAIR, "BCR-PDGFRA Fusion");

        assertNotNull(fusion);
        assertEquals("BCR", fusion.geneUp());
        assertEquals("PDGFRA", fusion.geneDown());
    }

    @Test
    public void ignoresFusionsOnUnknownGenes() {
        FusionExtractor fusionExtractor = testFusionExtractor();
        KnownFusionPair fusion = fusionExtractor.extract("IG", EventType.FUSION_PAIR, "IG-BCL2");

        assertNull(fusion);
    }

    @Test
    public void canExtractFusionPairsWithExonsUpDown() {
        FusionExtractor fusionExtractor = testFusionExtractor();
        KnownFusionPair fusion = fusionExtractor.extract("EGFR", EventType.FUSION_PAIR, "EGFRvII");

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
        FusionExtractor fusionExtractor =
                testFusionExtractorWithGeneChecker(new GeneChecker(Sets.newHashSet("IGH", "NKX2-1", "HLA-A", "ROS1")));
        KnownFusionPair fusion1 = fusionExtractor.extract("NKX2-1", EventType.FUSION_PAIR, "IGH-NKX2-1 Fusion");

        assertNotNull(fusion1);
        assertEquals("IGH", fusion1.geneUp());
        assertEquals("NKX2-1", fusion1.geneDown());

        KnownFusionPair fusion2 = fusionExtractor.extract("ROS1", EventType.FUSION_PAIR, "HLA-A-ROS1 Fusion");
        assertEquals("HLA-A", fusion2.geneUp());
        assertEquals("ROS1", fusion2.geneDown());

        KnownFusionPair fusion3 = fusionExtractor.extract("ROS1", EventType.FUSION_PAIR, "HLA-A-HLA-A");
        assertEquals("HLA-A", fusion3.geneUp());
        assertEquals("HLA-A", fusion3.geneDown());
    }

    @Test
    public void canExtractFusionPairsWithExons() {
        FusionExtractor fusionExtractor = testFusionExtractor();
        KnownFusionPair fusion = fusionExtractor.extract("MET", EventType.FUSION_PAIR_AND_EXON, "EXON 14 SKIPPING MUTATION");

        assertNotNull(fusion);
        assertEquals("MET", fusion.geneUp());
        assertEquals(13, (int) fusion.minExonUp());
        assertEquals(13, (int) fusion.maxExonUp());
        assertEquals("MET", fusion.geneDown());
        assertEquals(15, (int) fusion.minExonDown());
        assertEquals(15, (int) fusion.maxExonDown());
    }

    @Test
    public void canExtractExonicDelDupFusions() {
        FusionExtractor fusionExtractor = testFusionExtractorWithExonicDelDupKeyPhrases(Sets.newHashSet("skip this"));
        KnownFusionPair fusion = fusionExtractor.extract("EGFR", EventType.FUSION_PAIR, "skip this 20");

        assertNotNull(fusion);
        assertEquals("EGFR", fusion.geneUp());
        assertEquals(19, (int) fusion.minExonUp());
        assertEquals(19, (int) fusion.maxExonUp());
        assertEquals("EGFR", fusion.geneDown());
        assertEquals(21, (int) fusion.minExonDown());
        assertEquals(21, (int) fusion.maxExonDown());
    }

    @Test
    public void canFilterFusionPairsWithExonsOnWrongGenes() {
        FusionExtractor fusionExtractor = testFusionExtractor();
        assertNull(fusionExtractor.extract("BRAF", EventType.FUSION_PAIR_AND_EXON, "EXON 14 SKIPPING MUTATION"));
    }

    @Test
    public void canFilterNonConfiguredFusionPairsWithExons() {
        FusionExtractor fusionExtractor = testFusionExtractor();
        assertNull(fusionExtractor.extract("MET", EventType.FUSION_PAIR_AND_EXON, "Does not exist"));
    }

    @NotNull
    private static FusionExtractor testFusionExtractor() {
        return buildTestFusionExtractor(V37_GENE_CHECKER, Sets.newHashSet());
    }

    @NotNull
    private static FusionExtractor testFusionExtractorWithGeneChecker(@NotNull GeneChecker geneChecker) {
        return buildTestFusionExtractor(geneChecker, Sets.newHashSet());
    }

    @NotNull
    private static FusionExtractor testFusionExtractorWithExonicDelDupKeyPhrases(@NotNull Set<String> exonicDelDupKeyPhrases) {
        return buildTestFusionExtractor(V37_GENE_CHECKER, exonicDelDupKeyPhrases);
    }

    @NotNull
    private static FusionExtractor buildTestFusionExtractor(@NotNull GeneChecker geneChecker, @NotNull Set<String> exonicDelDupKeyPhrases) {
        return new FusionExtractor(geneChecker, new KnownFusionCache(), exonicDelDupKeyPhrases, true);
    }
}