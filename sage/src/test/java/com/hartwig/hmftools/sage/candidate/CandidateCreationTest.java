package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.sage.common.IndexedBases;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class CandidateCreationTest
{

    @Test
    public void testSoftClipInsert()
    {
        String refBaseStr = generateRandomBases(100);

        IndexedBases refBases = new IndexedBases(100, 0, refBaseStr.getBytes());

        String insertBases = "AAAAA";

        // first an insert on the right
        int scStartRefIndex = 51;
        String scBases = insertBases + refBaseStr.substring(scStartRefIndex, 71);
        String readBases = refBaseStr.substring(20, scStartRefIndex) + scBases;
        int scReadIndex = 31;

        AltRead altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, false);

        assertNotNull(altRead);
        assertEquals("T", altRead.Ref);
        assertEquals("T" + insertBases, altRead.Alt);

        // then on the left
        scStartRefIndex = 0;
        scBases = refBaseStr.substring(scStartRefIndex, 20) + insertBases;
        readBases = scBases + refBaseStr.substring(20, 51);
        scReadIndex = 0;

        altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, true);

        assertNotNull(altRead);
        assertEquals("C", altRead.Ref);
        assertEquals("C" + insertBases, altRead.Alt);

    }


}
