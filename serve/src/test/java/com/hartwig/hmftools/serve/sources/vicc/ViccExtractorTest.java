package com.hartwig.hmftools.serve.sources.vicc;

import org.junit.Ignore;
import org.junit.Test;

public class ViccExtractorTest {

    @Test
    @Ignore
    public void canExtractAnnotationForEntryWithTwoCodonFeatures() {
        // TODO Implement
//        CodonExtractor extractor = createWithDriverGenes(createDriverGenes("TP53", "PIK3CA", "KRAS"));
//
//        Feature feature1 = ViccTestFactory.testFeatureWithGeneAndName("PIK3CA", "E545X");
//        Feature feature2 = ViccTestFactory.testFeatureWithGeneAndName("KRAS", "G12X");
//        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(feature1, feature2));
//
//        assertEquals(2, extractor.extract(entry).size());
//        assertEquals(1, extractor.extract(entry).get(feature1).size());
//        assertEquals(1, extractor.extract(entry).get(feature2).size());
    }

    @Test
    @Ignore
    public void canExtractAnnotationForEntryWithTwoFeaturesExon() {
        // TODO Implement
//        ExonExtractor extractor = createWithDriverGenes(Lists.newArrayList());
//
//        Feature feature1 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions/deletions");
//        Feature feature2 = ViccTestFactory.testFeatureWithGeneAndName("ERBB2", "Exon 20 insertions");
//        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(feature1, feature2));
//
//        assertEquals(2, extractor.extract(entry).size());
//        assertEquals(1, extractor.extract(entry).get(feature1).size());
//        assertEquals(1, extractor.extract(entry).get(feature2).size());
    }

}