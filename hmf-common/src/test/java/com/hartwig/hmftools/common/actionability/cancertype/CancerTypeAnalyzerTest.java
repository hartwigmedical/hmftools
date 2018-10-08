package com.hartwig.hmftools.common.actionability.cancertype;

import static org.junit.Assert.assertFalse;

import com.google.common.collect.Lists;

import org.junit.Test;

public class CancerTypeAnalyzerTest {

    @Test
    public void foundDoids() {

        CancerTypeReading cancerTypeReading = ImmutableCancerTypeReading.builder()
                .cancerType("Lung")
                .doidSet("1324")
                .build();

        CancerTypeAnalyzer cancerTypeAnalyzer = new CancerTypeAnalyzer(Lists.newArrayList(cancerTypeReading));
        assertFalse(cancerTypeAnalyzer.foundTumorLocation("Skin Melenoma", "1324"));
   //     assertTrue(cancerTypeAnalyzer.foundTumorLocation("Lung Non-Small cell", "1324"));
    }
}