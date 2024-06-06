package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.sage.evidence.RawContext.exceedsSoftClipLowBaseQual;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class RefContextConsumerTest
{
    @Test
    public void testMnvLength()
    {
        assertEquals(1, RefContextConsumer.mnvLength(0, 0, "A".getBytes(), "CC".getBytes()));
        assertEquals(1, RefContextConsumer.mnvLength(0, 0, "AA".getBytes(), "C".getBytes()));

        assertEquals(2, RefContextConsumer.mnvLength(0, 0, "AAA".getBytes(), "CC".getBytes()));
        assertEquals(2, RefContextConsumer.mnvLength(0, 0, "AA".getBytes(), "CCC".getBytes()));
        assertEquals(2, RefContextConsumer.mnvLength(0, 0, "AAC".getBytes(), "CCC".getBytes()));

        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "AAA".getBytes(), "CCC".getBytes()));
        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "ACA".getBytes(), "CCCC".getBytes()));
        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "ACAC".getBytes(), "CCCC".getBytes()));

        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "AAAA".getBytes(), "CCCC".getBytes()));
        assertEquals(3, RefContextConsumer.mnvLength(1, 0, "AAAA".getBytes(), "CCCC".getBytes()));
        assertEquals(2, RefContextConsumer.mnvLength(2, 0, "AAAA".getBytes(), "CCCC".getBytes()));
        assertEquals(1, RefContextConsumer.mnvLength(3, 0, "AAAA".getBytes(), "CCCC".getBytes()));
    }

    @Test
    public void testLowQualSoftClips()
    {
        final byte[] baseQualities = new byte[30];

        for(int i = 0; i < baseQualities.length; ++i)
        {
            baseQualities[i] = (byte)MIN_SOFT_CLIP_MIN_BASE_QUAL;
        }

        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 0, 10));
        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 20, 10));

        for(int i = 5; i < 25; ++i)
        {
            baseQualities[i] = (byte)(MIN_SOFT_CLIP_MIN_BASE_QUAL - 1);
        }

        assertTrue(exceedsSoftClipLowBaseQual(baseQualities, 0, 10));
        assertTrue(exceedsSoftClipLowBaseQual(baseQualities, 20, 10));

        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 0, 7));
        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 23, 7));
    }
}
