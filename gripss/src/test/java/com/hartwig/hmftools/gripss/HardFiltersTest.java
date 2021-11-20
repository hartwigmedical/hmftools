package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createBreakend;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class HardFiltersTest
{
    private GripssTestApplication mGripss;

    public HardFiltersTest()
    {
        mGripss = new GripssTestApplication();
    }

    @Test
    public void testHardFilters()
    {
        VariantContext var1 = createBreakend(mGripss.IdGen.next(true), CHR_1, 100, "A", "G");

        mGripss.processVariant(var1);
    }

}
