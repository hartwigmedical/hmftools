package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MixedGermlineTest {

    @Test
    public void testNonCodingMnv() {
        testMixed(100, false);
    }

    @Test
    public void testCodingMnv() {
        testMixed(115251156, true);
        testMixed(115251157, false);
        testMixed(115251158, false);
        testMixed(115251159, true);
    }


    private void testMixed(int germlinePosition, boolean mvnPass) {
        final List<SageVariant> consumer = Lists.newArrayList();

        final SageVariant germlineSnv = createGermline(germlinePosition, "A", "G");
        final SageVariant somaticSnv = create(germlinePosition + 2, "A", "G");
        final SageVariant mixedMnv = create(germlinePosition, "ACA", "GCG");

        somaticSnv.localPhaseSet(3);
        mixedMnv.localPhaseSet(3);

        final Phase victim = new Phase(SageConfigTest.testConfig(), "1", consumer::add);

        victim.accept(mixedMnv);
        victim.accept(germlineSnv);
        victim.accept(somaticSnv);

        assertEquals(0, consumer.size());

        victim.flush();
        assertEquals(3, consumer.size());

        assertEquals(1, mixedMnv.mixedGermlineImpact());
        assertEquals(1, somaticSnv.mixedGermlineImpact());
        assertEquals(1, germlineSnv.mixedGermlineImpact());

        assertFalse(germlineSnv.isPassing());
        assertEquals(!mvnPass, somaticSnv.isPassing());
        assertEquals(mvnPass, mixedMnv.isPassing());
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
        ReadContextCounter counter = dummyCounter(variant, Strings.EMPTY);
        return new SageVariant(SageVariantTier.PANEL, variant, Sets.newHashSet(), Lists.newArrayList(), Lists.newArrayList(counter));
    }

    @NotNull
    static ReadContextCounter dummyCounter(@NotNull VariantHotspot variant, @NotNull final String microhomology) {
        ReadContext dummyReadContext = ReadContextCounterTest.readContext(100, 0, 0, 0, "AAA", microhomology);
        return new ReadContextCounter("SAMPLE",
                variant,
                dummyReadContext,
                new QualityRecalibrationMap(Collections.emptyList()),
                SageVariantTier.PANEL,
                1000,
                false);
    }

}
