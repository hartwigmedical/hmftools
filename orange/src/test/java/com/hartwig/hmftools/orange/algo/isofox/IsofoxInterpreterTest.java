package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.RnaFusion;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class IsofoxInterpreterTest {

    @Test
    public void canExtractGeneUpDownFromRNAFusion() {
        RnaFusion proper = IsofoxTestFactory.rnaFusionBuilder().name("X_Y").build();
        assertEquals("X", IsofoxInterpreter.geneUp(proper));
        assertEquals("Y", IsofoxInterpreter.geneDown(proper));

        RnaFusion upOnly = IsofoxTestFactory.rnaFusionBuilder().name("X_").build();
        assertEquals("X", IsofoxInterpreter.geneUp(upOnly));
        assertNull(IsofoxInterpreter.geneDown(upOnly));

        RnaFusion downOnly = IsofoxTestFactory.rnaFusionBuilder().name("_Y").build();
        assertNull(IsofoxInterpreter.geneUp(downOnly));
        assertEquals("Y", IsofoxInterpreter.geneDown(downOnly));

        RnaFusion none = IsofoxTestFactory.rnaFusionBuilder().name("_").build();
        assertNull(IsofoxInterpreter.geneUp(none));
        assertNull(IsofoxInterpreter.geneDown(none));

        RnaFusion empty = IsofoxTestFactory.rnaFusionBuilder().name(Strings.EMPTY).build();
        assertNull(IsofoxInterpreter.geneUp(empty));
        assertNull(IsofoxInterpreter.geneDown(empty));
    }
}