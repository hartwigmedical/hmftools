package com.hartwig.hmftools.purple.somatic;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class KataegisWindowTest
{
    @Test
    public void testAvgDistance()
    {
        final SomaticVariant context1 = KataegisQueueTest.create("1", 1000, false);
        final SomaticVariant context2 = KataegisQueueTest.create("1", 1100, false);

        final KataegisWindow window = new KataegisWindow(context1);
        assertEquals(0, window.count());

        window.add(context1);
        assertEquals(1, window.count());
        assertEquals(0, window.averageDistance());

        window.add(context2);
        assertEquals(2, window.count());
        assertEquals(100, window.averageDistance());
    }

    @Test
    public void testAvgDistanceRounding()
    {
        final KataegisWindow window = new KataegisWindow(KataegisQueueTest.create("1", 54730299, false));
        for(int i = 0; i < 177; i++)
        {
            window.add(KataegisQueueTest.create("1", 54906422, false));
        }
        assertEquals(1001, window.averageDistance());
    }
}
