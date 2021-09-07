package com.hartwig.hmftools.purple.segment;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.junit.Test;

public class StructuralVariantPositionFactoryTest
{

    @Test
    public void excludeInserts()
    {
        final StructuralVariant variant =
                PurpleTestUtils.createStructuralVariant("1", 1001, "1", 1001, StructuralVariantType.INS).build();

        assertEquals(0, SVSegmentFactory.create(Lists.newArrayList(variant)).size());
    }
}
