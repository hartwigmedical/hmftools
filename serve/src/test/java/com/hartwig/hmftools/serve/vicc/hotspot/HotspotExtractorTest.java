package com.hartwig.hmftools.serve.vicc.hotspot;

import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotExtractorTest {

    @Test
    public void canAssessWhetherFeatureIsProteinAnnotation() {
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("L2230V"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("V5del"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("L755_T759del"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("G10dup"));

        assertFalse(HotspotExtractor.isResolvableProteinAnnotation(Strings.EMPTY));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("truncating"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("20LtoV"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L20"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("LP"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L2"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L20Pdel5"));

        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("963_D1010splice"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("T599insTT"));
    }
}