package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneDown;
import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneUp;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.rna.RnaFusion;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class FusionNameUtilTest
{
    @Test
    public void canExtractGeneUpDownFromRnaFusion()
    {
        RnaFusion proper = IsofoxTestFactory.rnaFusionBuilder().name("X_Y").build();
        assertEquals("X", geneUp(proper));
        assertEquals("Y", geneDown(proper));

        RnaFusion upOnly = IsofoxTestFactory.rnaFusionBuilder().name("X_").build();
        assertEquals("X", geneUp(upOnly));
        assertNull(geneDown(upOnly));

        RnaFusion downOnly = IsofoxTestFactory.rnaFusionBuilder().name("_Y").build();
        assertNull(geneUp(downOnly));
        assertEquals("Y", geneDown(downOnly));

        RnaFusion none = IsofoxTestFactory.rnaFusionBuilder().name("_").build();
        assertNull(geneUp(none));
        assertNull(geneDown(none));

        RnaFusion empty = IsofoxTestFactory.rnaFusionBuilder().name(Strings.EMPTY).build();
        assertNull(geneUp(empty));
        assertNull(geneDown(empty));
    }
}