package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createStructuralVariant;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;
import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.junit.Test;

public class StructuralVariantPositionFactoryTest
{
    @Test
    public void excludeInserts()
    {
        StructuralVariant variant = createStructuralVariant("1", 1001, "1", 1001, INS).build();
        assertEquals(0, SVSegmentFactory.create(Lists.newArrayList(variant)).size());
    }

    @Test
    public void deletions()
    {
        StructuralVariant variant = createStructuralVariant("1", 1001, "1", 2001, DEL).build();
        List<SVSegment> segments = SVSegmentFactory.create(of(variant));
        assertEquals(2, segments.size());
        assertEquals(1002, segments.get(0).position());
        assertEquals(2001, segments.get(1).position());
    }

    @Test
    public void inversions()
    {
        StructuralVariant variant = createStructuralVariant("1", 1001, "1", 2001, INV).build();
        List<SVSegment> segments = SVSegmentFactory.create(of(variant));
        assertEquals(2, segments.size());
        assertEquals(1002, segments.get(0).position());
        assertEquals(2002, segments.get(1).position());
    }
}