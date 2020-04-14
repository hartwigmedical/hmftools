package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfigTest;
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextCounterTest;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MixedGermlineTest {

    @Test
    public void testStandard() {
        final List<SageVariant> consumer = Lists.newArrayList();

        final SageVariant germlineSnv = createGermline(100, "A", "G");
        final SageVariant somaticSnv = create(102, "A", "G");
        final SageVariant mixedMnv = create(100, "ACA", "GAG");

        somaticSnv.localPhaseSet(3);
        mixedMnv.localPhaseSet(3);

        final Phase victim = new Phase(SageConfigTest.testConfig(), consumer::add);

        victim.accept(mixedMnv);
        victim.accept(germlineSnv);
        victim.accept(somaticSnv);

        assertEquals(0, consumer.size());

        victim.flush();
        assertEquals(3, consumer.size());

        assertTrue(mixedMnv.filters().contains(SoftFilter.MIXED_GERMLINE_SOMATIC_MNV.toString()));
        assertEquals(1, mixedMnv.mixedGermlineEffect());
        assertEquals(1, somaticSnv.mixedGermlineEffect());
        assertEquals(0, germlineSnv.mixedGermlineEffect());
    }

    @NotNull
    private static SageVariant createGermline(long position, @NotNull String ref, @NotNull String alt) {
        SageVariant result = create(position, ref, alt);
        result.filters().add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.toString());
        return result;
    }

    @NotNull
    private static SageVariant create(long position, @NotNull String ref, @NotNull String alt) {
        VariantHotspot variant = LocalPhaseSetTest.create(position, ref, alt);
        ReadContext dummyReadContext = ReadContextCounterTest.readContext(100, 0, 0, 0, "AAA");
        ReadContextCounter counter = new ReadContextCounter("SAMPLE",
                variant,
                dummyReadContext,
                new QualityRecalibrationMap(Collections.emptyList()),
                SageVariantTier.PANEL,
                1000,
                false);

        return new SageVariant(SageVariantTier.PANEL, variant, Sets.newHashSet(), Lists.newArrayList(), Lists.newArrayList(counter));
    }

}
