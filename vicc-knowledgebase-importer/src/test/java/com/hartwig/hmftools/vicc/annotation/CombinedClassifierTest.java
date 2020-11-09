package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CombinedClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsCombinedEvent() {
        assertTrue(CombinedClassifier.isCombinedEvent("Ex19 del L858R", "EGFR"));

        assertTrue(CombinedClassifier.isCombinedEvent("JAK1 (S646F;R683)", "JAK1"));
        assertTrue(CombinedClassifier.isCombinedEvent("KIT (627-664,664-714,449-514)", "KIT"));
        assertFalse(CombinedClassifier.isCombinedEvent("KIT mutation in exon 9,11,13,14 or 17", "KIT"));

        assertTrue(CombinedClassifier.isCombinedEvent("EGFR inframe deletion (L747),inframe insertion (P753PS)", "EGFR"));
        assertTrue(CombinedClassifier.isCombinedEvent("KIT inframe deletion (416-422),inframe insertion (416-422)", "KIT"));
        assertTrue(CombinedClassifier.isCombinedEvent("CSF3R frameshift variant (D771),frameshift variant (S783)", "CSF3R"));
        assertTrue(CombinedClassifier.isCombinedEvent("NOTCH1 splice donor variant (2245-2536),splice acceptor variant (2245-2536),"
                + "stop gained (2245-2536),stop lost (2245-2536),frameshift variant (2245-2536)", "NOTCH1"));

        assertTrue(CombinedClassifier.isCombinedEvent("BCR-ABL F486S", "ABL"));
        assertTrue(CombinedClassifier.isCombinedEvent("NPM1-ALK  amp", "ALK"));

        assertFalse(CombinedClassifier.isCombinedEvent("Exon 20 insertions/deletions", "ERBB2"));
        assertFalse(CombinedClassifier.isCombinedEvent("3'UTR alteration (c.642+70C>A)", "VHL"));
    }
}