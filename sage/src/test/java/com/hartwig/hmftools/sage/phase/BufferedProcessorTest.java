package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.junit.Test;

public class BufferedProcessorTest {

    @Test
    public void testLongerContainsShorter() {
        final VariantHotspot mnv = create(100, "CAC", "TGT");

        assertTrue(BufferedProcessor.longerContainsShorter(create(100, "C", "T"), mnv));
        assertFalse(BufferedProcessor.longerContainsShorter(create(100, "C", "C"), mnv));

        assertTrue(BufferedProcessor.longerContainsShorter(create(101, "A", "G"), mnv));
        assertTrue(BufferedProcessor.longerContainsShorter(create(102, "C", "T"), mnv));

        assertTrue(BufferedProcessor.longerContainsShorter(create(100, "CA", "TG"), mnv));
        assertTrue(BufferedProcessor.longerContainsShorter(create(101, "AC", "GT"), mnv));

        assertTrue(BufferedProcessor.longerContainsShorter(create(100, "CAC", "TGT"), mnv));
    }

    static VariantHotspot create(long position, String ref, String alt) {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").position(position).ref(ref).alt(alt).build();
    }

}
