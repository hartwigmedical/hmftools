package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.serve.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Test;

public class FusionExtractorTest {

    @Test
    public void canExtractFusionPairsGenesUnknown() {
        FusionExtractor fusionExtractor = new FusionExtractor(HmfGenePanelSupplier.allGenesMap37());
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "IG",
                "IG-BCL2",
                "description",
                "chromosome",
                "pos", null);
        assertEquals(fusionsPerFeature, fusionExtractor.extractFusionPairs(viccEntry));
    }

    @Test
    public void canExtractFusionPairsGenes() {
        FusionExtractor fusionExtractor = new FusionExtractor(HmfGenePanelSupplier.allGenesMap37());
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "PDGFRA",
                "BCR-PDGFRA Fusion",
                "description",
                "chromosome",
                "pos", null);
        fusionsPerFeature.put(viccEntry.features().get(0), ImmutableKnownFusionPair.builder().geneUp("BCR").geneDown("PDGFRA").build());
        assertEquals(fusionsPerFeature, fusionExtractor.extractFusionPairs(viccEntry));
    }

    @Test
    public void canExtractFusionPairsWithExonsUpDown() {
        FusionExtractor fusionExtractor = new FusionExtractor(HmfGenePanelSupplier.allGenesMap37());
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "EGFR",
                "EGFRvII",
                "description",
                "chromosome",
                "pos", null);
        fusionsPerFeature.put(viccEntry.features().get(0),
                ImmutableKnownFusionPair.builder()
                        .geneUp("EGFR")
                        .maxExonUp(13)
                        .minExonUp(13)
                        .geneDown("EGFR")
                        .maxExonDown(16)
                        .minExonDown(16)
                        .build());
        assertEquals(fusionsPerFeature, fusionExtractor.extractFusionPairs(viccEntry));
    }

    @Test
    public void canExtractFusionPairsWithRangeExons() {
        FusionExtractor fusionExtractor = new FusionExtractor(HmfGenePanelSupplier.allGenesMap37());
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC,
                "any",
                "MET",
                "EXON 14 SKIPPING MUTATION",
                "MET EXON 14 SKIPPING MUTATION",
                "chromosome",
                "pos", null);
        fusionsPerFeature.put(viccEntry.features().get(0),
                ImmutableKnownFusionPair.builder()
                        .geneUp("MET")
                        .maxExonUp(13)
                        .minExonUp(13)
                        .geneDown("MET")
                        .maxExonDown(15)
                        .minExonDown(15)
                        .build());
        fusionExtractor.extractFusionPairs(viccEntry);
        assertEquals(fusionsPerFeature, fusionExtractor.extractFusionPairs(viccEntry));
    }

}