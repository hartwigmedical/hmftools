package com.hartwig.hmftools.common.actionability.util;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MultiDrugCuratorTest {

    @Test
    public void canReformatDrugs() {
        assertEquals("drugA", MultiDrugCurator.reformat("drugA"));
        assertEquals("drugA + drugB", MultiDrugCurator.reformat("drugA + drugB"));
        assertEquals("drugA + drugB", MultiDrugCurator.reformat("drugB + drugA"));
    }
}