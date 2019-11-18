package com.hartwig.hmftools.strelka.mnv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GapReadsTest {

    @Test
    public void testFavourActualBases() {
        final GapReads empty = GapReads.empty();
        assertEquals((Character) 'N', empty.mostFrequentRead());

        final GapReads multiN =  GapReads.addRead(GapReads.addRead(empty, 'N'), 'N');
        assertEquals((Character) 'N', multiN.mostFrequentRead());

        final GapReads multiNSingleA = GapReads.addRead(multiN, 'A');
        assertEquals((Character) 'A', multiNSingleA.mostFrequentRead());
    }
}
