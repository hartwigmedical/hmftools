package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.Microhomology.findHomology;
import static com.hartwig.hmftools.sage.common.Microhomology.findLeftHomologyShift;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class HomologyTest
{
    @Test
    public void testHomologyLeftAlignment()
    {
        //            60        70
        //            01234567890123
        // ref bases: CTCGGAAAAAAAAC

        // insert scenario:
        SimpleVariant var = createSimpleVariant(70, "A", "AAAAA");

        String readBases = REF_BASES_200.substring(60, 70) + "AAAAA" + REF_BASES_200.substring(71, 80);
        int leftHomShift = findLeftHomologyShift(var, REF_SEQUENCE_200, readBases.getBytes(), 10);
        assertEquals(6, leftHomShift);

        // del scenario
        var = createSimpleVariant(68, "AAAAA", "A");

        readBases = REF_BASES_200.substring(60, 69) + REF_BASES_200.substring(73, 80);
        leftHomShift = findLeftHomologyShift(var, REF_SEQUENCE_200, readBases.getBytes(), 8);
        assertEquals(4, leftHomShift);

    }

    @Test
    public void testHomology()
    {
        SimpleVariant var = createSimpleVariant(26, "A", "AT");

        String readBases = REF_BASES_200.substring(10, 26) + "AT" + REF_BASES_200.substring(27, 40);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "16M1I13M";
        SAMRecord read = buildSamRecord(10, readCigar, readBases, baseQuals);

        Microhomology homology = findHomology(var, read.getReadBases(), 16);

        assertNotNull(homology);
        assertEquals("T", homology.Bases);
        assertEquals(4, homology.Length);

        var = createSimpleVariant(26, "A", "ATT");

        readBases = REF_BASES_200.substring(10, 26) + "ATT" + REF_BASES_200.substring(27, 40);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "16M2I13M";
        read = buildSamRecord(10, readCigar, readBases, baseQuals);

        homology = findHomology(var, read.getReadBases(), 16);

        assertNotNull(homology);
        assertEquals("TT", homology.Bases);
        assertEquals(4, homology.Length);

        // delete
        var = createSimpleVariant(64, "GAAA", "G");

        readBases = REF_BASES_200.substring(50, 65) + REF_BASES_200.substring(69, 80);
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "15M3D11M";
        read = buildSamRecord(50, readCigar, readBases, baseQuals);

        homology = findHomology(var, read.getReadBases(), 15);

        assertNotNull(homology);
        assertEquals("AAA", homology.Bases);
        assertEquals(3, homology.Length);

        // checks read bases not ref bases to determine homology
        var = createSimpleVariant(26, "A", "AAAA");
        readBases = REF_BASES_200.substring(10, 26) + "AAAAAAAAAAAAAAAAAAAAATG";
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "16M3I20M";
        read = buildSamRecord(10, readCigar, readBases, baseQuals);

        homology = findHomology(var, read.getReadBases(), 16);

        assertNotNull(homology);
        assertEquals("AAA", homology.Bases);
        assertEquals(17, homology.Length);

        // multiple copies and then finishes with a partial copy
        // eg G(AACTC)AACTCAACTCAACCCTTT -> GAACTCAACTCAA(CTCAA)CCCTTT

        var = createSimpleVariant(26, "G", "GAACTC");
        readBases = REF_BASES_200.substring(10, 26) + "GAACTCAACTCAACTCAACCCTTT";
        baseQuals = buildDefaultBaseQuals(readBases.length());
        readCigar = "16M5I18M";
        read = buildSamRecord(10, readCigar, readBases, baseQuals);

        homology = findHomology(var, read.getReadBases(), 16);

        assertNotNull(homology);
        assertEquals("AACTC", homology.Bases);
        assertEquals(13, homology.Length);
    }
}
