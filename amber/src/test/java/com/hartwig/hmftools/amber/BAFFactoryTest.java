package com.hartwig.hmftools.amber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.pileup.PileupFile;

import org.junit.Test;

public class BAFFactoryTest {

    @Test
    public void testDepths() {
        final Pileup normal = PileupFile.fromString("seq2\t156\tA\t11\t.$......+2AG.+2AG.+2AGGG\t<975;:<<<<<");
        final Pileup tumor = PileupFile.fromString("seq2\t156\tA\t31\t.$......+2AG.+2AG.+2AGGG\t<975;:<<<<<");
        assertNotEquals(normal.readCount(), tumor.readCount());

        final AmberBAF victim = BAFFactory.create('A', normal, tumor);
        assertEquals(normal.readCount(), victim.normalDepth());
        assertEquals(tumor.readCount(), victim.tumorDepth());
    }
}
