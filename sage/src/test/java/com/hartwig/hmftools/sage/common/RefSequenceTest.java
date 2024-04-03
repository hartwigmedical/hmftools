package com.hartwig.hmftools.sage.common;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class RefSequenceTest
{
    @Test
    public void testRefSequence()
    {
        String refBases = MockRefGenome.generateRandomBases(101);
        RefSequence refSequence = new RefSequence(100, refBases.getBytes());

        assertEquals(200, refSequence.End);
        assertEquals(10, refSequence.index(110));
        assertEquals(101, refSequence.length());
        assertEquals(refBases.substring(5, 16), refSequence.positionBases(105, 115));
        assertEquals(refBases.substring(5, 16), refSequence.indexBases(5, 15));
        assertEquals(refBases.substring(9, 12), new String(refSequence.trinucleotideContext(110)));
    }
}
