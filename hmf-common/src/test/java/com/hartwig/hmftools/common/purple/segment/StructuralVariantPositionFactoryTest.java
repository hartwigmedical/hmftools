package com.hartwig.hmftools.common.purple.segment;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.junit.Test;

public class StructuralVariantPositionFactoryTest {

    @Test
    public void excludeInserts() {
        final StructuralVariant variant =
                PurpleDatamodelTest.createStructuralVariant("1", 1001, "1", 1001, StructuralVariantType.INS).build();

        assertEquals(0, SVSegmentFactory.create(Lists.newArrayList(variant)).size());
    }
}
