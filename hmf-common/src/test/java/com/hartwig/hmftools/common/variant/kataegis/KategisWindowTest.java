package com.hartwig.hmftools.common.variant.kataegis;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class KategisWindowTest {

    @Test
    public void testAvgDistance() {
        final VariantContext context1 = KataegisQueueTest.create("1", 1000, false);
        final VariantContext context2 = KataegisQueueTest.create("1", 1100, false);

        final KataegisWindow window = new KataegisWindow(context1);
        assertEquals(0, window.count());

        window.add(context1);
        assertEquals(1, window.count());
        assertEquals(0, window.averageDistance());

        window.add(context2);
        assertEquals(2, window.count());
        assertEquals(100, window.averageDistance());

    }

}
