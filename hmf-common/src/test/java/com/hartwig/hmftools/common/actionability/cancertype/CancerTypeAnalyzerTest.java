package com.hartwig.hmftools.common.actionability.cancertype;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.junit.Test;

public class CancerTypeAnalyzerTest {

    @Test
    public void canFindDoids() {
        CancerTypeReading cancerTypeReading = ImmutableCancerTypeReading.builder()
                .cancerType("Lung")
                .doidSet("1324")
                .build();

        CancerTypeAnalyzer cancerTypeAnalyzer = new CancerTypeAnalyzer(Lists.newArrayList(cancerTypeReading));
        assertFalse(cancerTypeAnalyzer.foundTumorLocation("Skin Melanoma", "1324"));
        assertTrue(cancerTypeAnalyzer.foundTumorLocation("Lung", "1324"));
        assertFalse(cancerTypeAnalyzer.foundTumorLocation("Lung Cancer: Non-Small Cell", "1324"));
    }
}