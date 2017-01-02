package com.hartwig.hmftools.patientreporter.slicing;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SlicerTest {

    @Test
    public void includeCheckWorksAsExpected() {
        final Slicer slicer = buildSlicer();
        final SomaticVariant include1 = buildVariant("X", 150);
        final SomaticVariant include2 = buildVariant("X", 300);
        final SomaticVariant include3 = buildVariant("X", 400);
        final SomaticVariant include4 = buildVariant("Y", 570);

        final SomaticVariant notInclude1 = buildVariant("X", 1);
        final SomaticVariant notInclude2 = buildVariant("X", 250);
        final SomaticVariant notInclude3 = buildVariant("Y", 1000);
        final SomaticVariant notInclude4 = buildVariant("NotExists", 100);

        assertTrue(slicer.includes(include1));
        assertTrue(slicer.includes(include2));
        assertTrue(slicer.includes(include3));
        assertTrue(slicer.includes(include4));

        assertFalse(slicer.includes(notInclude1));
        assertFalse(slicer.includes(notInclude2));
        assertFalse(slicer.includes(notInclude3));
        assertFalse(slicer.includes(notInclude4));
    }

    @NotNull
    private static Slicer buildSlicer() {
        SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put("X", new GenomeRegion(100, 200));
        regionMap.put("X", new GenomeRegion(300, 400));
        regionMap.put("Y", new GenomeRegion(500, 600));
        return new Slicer(regionMap);
    }

    @NotNull
    private static SomaticVariant buildVariant(@NotNull String chromosome, long position) {
        return new SomaticVariant(VariantType.SNP, "Some Filter", chromosome, position, Lists.newArrayList(), 0D,
                false, false);
    }
}
