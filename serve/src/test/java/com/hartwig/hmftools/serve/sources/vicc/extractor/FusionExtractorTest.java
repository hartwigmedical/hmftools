package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.junit.Test;

public class FusionExtractorTest {

    @Test
    public void canExtractFusionPairsGenes() {
        FusionExtractor fusionExtractor = new FusionExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("PDGFRA", "BCR-PDGFRA Fusion");

        Map<Feature, KnownFusionPair> fusionsPerFeature = fusionExtractor.extractFusionPairs(viccEntry);
        assertEquals(1, fusionsPerFeature.size());
        assertEquals("BCR", fusionsPerFeature.get(viccEntry.features().get(0)).geneUp());
        assertEquals("PDGFRA", fusionsPerFeature.get(viccEntry.features().get(0)).geneDown());
    }

    @Test
    public void ignoresFusionsOnUnknownGenes() {
        FusionExtractor fusionExtractor = new FusionExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("IG", "IG-BCL2");

        Map<Feature, KnownFusionPair> fusionsPerFeature = fusionExtractor.extractFusionPairs(viccEntry);

        assertTrue(fusionsPerFeature.isEmpty());
    }

    @Test
    public void canExtractFusionPairsWithExonsUpDown() {
        FusionExtractor fusionExtractor = new FusionExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("EGFR", "EGFRvII");

        Map<Feature, KnownFusionPair> fusionsPerFeature = fusionExtractor.extractFusionPairs(viccEntry);
        assertEquals(1, fusionsPerFeature.size());
        assertEquals("EGFR", fusionsPerFeature.get(viccEntry.features().get(0)).geneUp());
        assertEquals(13, (int) fusionsPerFeature.get(viccEntry.features().get(0)).minExonUp());
        assertEquals(13, (int) fusionsPerFeature.get(viccEntry.features().get(0)).maxExonUp());
        assertEquals("EGFR", fusionsPerFeature.get(viccEntry.features().get(0)).geneDown());
        assertEquals(16, (int) fusionsPerFeature.get(viccEntry.features().get(0)).minExonDown());
        assertEquals(16, (int) fusionsPerFeature.get(viccEntry.features().get(0)).maxExonDown());
    }

    @Test
    public void canExtractFusionPairsWithOddNames() {
        FusionExtractor fusionExtractor = new FusionExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("NKX2-1", "IGH-NKX2-1 Fusion");

        Map<Feature, KnownFusionPair> fusionsPerFeature = fusionExtractor.extractFusionPairs(viccEntry);
        assertEquals(1, fusionsPerFeature.size());
        assertEquals("IGH", fusionsPerFeature.get(viccEntry.features().get(0)).geneUp());
        assertEquals("NKX2-1", fusionsPerFeature.get(viccEntry.features().get(0)).geneDown());
    }

    @Test
    public void canExtractFusionPairsWithRangeExons() {
        FusionExtractor fusionExtractor = new FusionExtractor(GeneChecker.buildForHG19());
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("MET", "EXON 14 SKIPPING MUTATION");

        Map<Feature, KnownFusionPair> fusionsPerFeature = fusionExtractor.extractFusionPairs(viccEntry);
        assertEquals(1, fusionsPerFeature.size());
        assertEquals("MET", fusionsPerFeature.get(viccEntry.features().get(0)).geneUp());
        assertEquals(13, (int) fusionsPerFeature.get(viccEntry.features().get(0)).minExonUp());
        assertEquals(13, (int) fusionsPerFeature.get(viccEntry.features().get(0)).maxExonUp());
        assertEquals("MET", fusionsPerFeature.get(viccEntry.features().get(0)).geneDown());
        assertEquals(15, (int) fusionsPerFeature.get(viccEntry.features().get(0)).minExonDown());
        assertEquals(15, (int) fusionsPerFeature.get(viccEntry.features().get(0)).maxExonDown());
    }
}