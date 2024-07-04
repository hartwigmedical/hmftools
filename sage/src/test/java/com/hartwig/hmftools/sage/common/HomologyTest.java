package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.Microhomology.findHomology;
import static com.hartwig.hmftools.sage.common.Microhomology.findLeftHomologyShift;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.candidate.RefContextCache;
import com.hartwig.hmftools.sage.candidate.RefContextConsumer;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

public class HomologyTest
{
    @Test
    public void testHomologyLeftAlignment()
    {
        //            50        60        70        80
        //            0123456789012345678901234567890
        // ref bases: CTGTCTGTGACTCGGAAAAAAAACTCCCTGA

        // test 1: insert scenario:
        SimpleVariant var = createSimpleVariant(70, "A", "AAAAA");

        String readBases = REF_BASES_200.substring(60, 70) + "AAAAA" + REF_BASES_200.substring(71, 80);
        int leftHomShift = findLeftHomologyShift(var, REF_SEQUENCE_200, readBases.getBytes(), 10);
        assertEquals(6, leftHomShift);

        // test 2: dinucleotide repeats
        String refBases = "CTGTCTGTGACTCGGATATATATATATATATCCCCTTGCGCTTCCCAGGT";
        RefSequence refSequence = new RefSequence(0, refBases.getBytes());
        //            10        20        30
        //            012345678901234567890123
        // ref bases: CTCGGATATATATATATATATCCC

        var = createSimpleVariant(30, "T", "TATAT");

        readBases = refBases.substring(10, 30) + "TATAT" + refBases.substring(31, 40);
        leftHomShift = findLeftHomologyShift(var, refSequence, readBases.getBytes(), 20);
        assertEquals(16, leftHomShift);

        // test 3: duplicate scenario
        refBases = "CTGTCTGTGACAAACCCGGGTCGGATCCCGGTAGGTAT";
        //          0         10        20        30
        //          0123456789012345678901234567890123456789
        refSequence = new RefSequence(0, refBases.getBytes());

        String insertedBases = "AAACCCGGG";
        var = createSimpleVariant(19, "G", "G" + insertedBases);

        readBases = refBases.substring(10, 19) + "G" + insertedBases + refBases.substring(20, 35);

        //          0         10        20        30        40
        //          012345678901234567890123456789012345678901234567890
        // read     CTGTCTGTGACAAACCCGGGAAACCCGGGTCGGATCCCGGTAGGTAT

        leftHomShift = findLeftHomologyShift(var, refSequence, readBases.getBytes(), 9);
        assertEquals(9, leftHomShift);

        // test 4: a duplication later on
        refBases = "CTGTCTGTGACAAACCCGGGAAACCCGGGTCGGATCCCGGTAGGTAT";
        //          0         10        20        30        40
        //          01234567890123456789012345678901234567890123456789
        refSequence = new RefSequence(0, refBases.getBytes());
        var = createSimpleVariant(28, "G", "G" + insertedBases);

        readBases = refBases.substring(10, 28) + "G" + insertedBases + refBases.substring(20, 35);

        leftHomShift = findLeftHomologyShift(var, refSequence, readBases.getBytes(), 18);
        assertEquals(18, leftHomShift);
    }

    @Test
    public void testHomologyLeftAlignmentImperfectHomology()
    {
        //                          10        20        30        40        50        60        70
        //                 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        String refBases = "TTTGGTGAGGAATGGGATACTGTCTGAAAGTATCTCTACAAAGGGATAAAGAGTAACTTTACATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAATAGGACAAATCA";
        String readBases = "CATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAATAGGATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAAT";
        String alt = "ATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAATAGGA";

        SimpleVariant var = createSimpleVariant(124, alt.substring(0, 1), alt);

        RefSequence refSequence = new RefSequence(1, refBases.getBytes());

        int leftHomShift = findLeftHomologyShift(var, refSequence, readBases.getBytes(), 62);
        assertEquals(62, leftHomShift);
    }

    @Test
    public void testHomologyLeftAlignmentCandidate()
    {
        // test left-aligning a candidate variant found in a right soft-clip
        // the soft-clip must be long enough to have the duplicated section and 12-bases of matching ref

        String refBases = "CTGTCTGTGACACGTTGCAGTCGGATCCCGGTAGGTATTGTGACAAACCGTAGCTTGACCAG";
        //                 0         10        20        30        40        50        60
        //                 01234567890123456789012345678901234567890123456789012345678901
        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        // up to and including the ref base are duplicated
        int varPosition = 30;
        int dupLength = 10;
        String duplicatedBases = refBases.substring(varPosition - dupLength + 1, varPosition + 1);

        String softClipBases = duplicatedBases + refBases.substring(31, 46); // then continues in the ref sequence

        int readStart = 5;
        String refAlignedBases = refBases.substring(readStart, varPosition + 1);
        String readBases = refAlignedBases + softClipBases;

        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 0, 200);
        RefContextCache refContextCache = new RefContextCache(TEST_CONFIG, Collections.emptyList(), Collections.emptyList());
        RefContextConsumer refContextConsumer = new RefContextConsumer(TEST_CONFIG, region, refSequence, refContextCache, Collections.emptyList());

        String cigar = buildCigarString(refAlignedBases.length(), 0, softClipBases.length());

        SAMRecord read = buildSamRecord(readStart, cigar, readBases);
        refContextConsumer.processRead(read);
        refContextConsumer.processRead(read);

        List<AltContext> altContexts = refContextCache.altContexts();
        TestCase.assertEquals(1, altContexts.size());

        AltContext altContext = altContexts.get(0);
        VariantReadContext readContext = altContext.readContext();
        assertEquals(20, readContext.variant().Position);

        // again but this time with an indel in the read bases to confuse where the implied ref position is
        int delPosition = 20;
        readBases = refAlignedBases.substring(0, delPosition) + refAlignedBases.substring(delPosition + 2) + softClipBases;

        cigar = "20M2D4M25S";
        read = buildSamRecord(readStart, cigar, readBases);

        refContextCache = new RefContextCache(TEST_CONFIG, Collections.emptyList(), Collections.emptyList());
        refContextConsumer = new RefContextConsumer(TEST_CONFIG, region, refSequence, refContextCache, Collections.emptyList());
        refContextConsumer.processRead(read);
        refContextConsumer.processRead(read);

        altContexts = refContextCache.altContexts();
        TestCase.assertEquals(2, altContexts.size());

        altContext = altContexts.stream().filter(x -> x.position() == 18).findFirst().orElse(null);
        assertNotNull(altContext);
        readContext = altContext.readContext();
        assertEquals(17, readContext.CorePositionStart);
        assertEquals(32, readContext.CorePositionEnd);
    }

    private static Microhomology findHomology(final SimpleVariant variant, final byte[] readBases, int varReadIndex)
    {
        return Microhomology.findHomology(variant, readBases, varReadIndex, true);
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
