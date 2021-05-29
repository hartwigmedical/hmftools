package com.hartwig.hmftools.serve.extraction.hotspot.tools;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class AnnotatedHotspotVCFCheckerTest {

    @Test
    public void canSmartMatch() {
        // Hotspots should never smart-match
        assertFalse(AnnotatedHotspotVCFChecker.isApproximateIndelMatch("p.K618I", "p.K617I"));

        // This is accepted
       assertTrue (AnnotatedHotspotVCFChecker.isApproximateIndelMatch("p.K618del", "p.K617del"));
    }

}