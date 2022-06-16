package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.RnaFusion;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class RNAFusionSelectorTest {

    @Test
    public void canExtractGeneUpDownFromRNAFusion() {
        RnaFusion proper = IsofoxTestFactory.rnaFusionBuilder().name("X_Y").build();
        assertEquals("X", RNAFusionSelector.geneUp(proper));
        assertEquals("Y", RNAFusionSelector.geneDown(proper));

        RnaFusion upOnly = IsofoxTestFactory.rnaFusionBuilder().name("X_").build();
        assertEquals("X", RNAFusionSelector.geneUp(upOnly));
        assertNull(RNAFusionSelector.geneDown(upOnly));

        RnaFusion downOnly = IsofoxTestFactory.rnaFusionBuilder().name("_Y").build();
        assertNull(RNAFusionSelector.geneUp(downOnly));
        assertEquals("Y", RNAFusionSelector.geneDown(downOnly));

        RnaFusion none = IsofoxTestFactory.rnaFusionBuilder().name("_").build();
        assertNull(RNAFusionSelector.geneUp(none));
        assertNull(RNAFusionSelector.geneDown(none));

        RnaFusion empty = IsofoxTestFactory.rnaFusionBuilder().name(Strings.EMPTY).build();
        assertNull(RNAFusionSelector.geneUp(empty));
        assertNull(RNAFusionSelector.geneDown(empty));
    }
}