package com.hartwig.hmftools.sage.candidate;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.candidate.RefContextConsumer;

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
}
