package com.hartwig.hmftools.serve.vicc.hotspot;

import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotExtractorTest {

    @Test
    public void canAssessWhetherFeatureIsProteinAnnotation() {
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("K5N"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("L2230V"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("V5del"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("L755_T759del"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("L755_T756delinsPP"));
        assertTrue(HotspotExtractor.isResolvableProteinAnnotation("G10dup"));

        // Just plain wrong protein annotations
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation(Strings.EMPTY));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("truncating"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("20LtoV"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L20"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("LP"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L2"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L20Pdel5"));

        // Splice variants are ignored by hotspot extractor:
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("963_D1010splice"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("T599insTT"));

        // Inframe event is too long
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L4_T40del"));
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("L698_S1037dup"));

        // Wild-type mutations are ignored
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(HotspotExtractor.isResolvableProteinAnnotation("V600E/K"));
    }
}