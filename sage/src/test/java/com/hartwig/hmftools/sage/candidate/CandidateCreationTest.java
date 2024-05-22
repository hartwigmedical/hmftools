package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.TEST_LEFT_FLANK;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class CandidateCreationTest
{
    @Test
    public void testSoftClipInsert()
    {
        String refBaseStr = generateRandomBases(100);

        RefSequence refBases = new RefSequence(100, refBaseStr.getBytes());

        String insertBases = "AAAAA";

        // first an insert on the right
        int scStartRefIndex = 51;
        String scBases = insertBases + refBaseStr.substring(scStartRefIndex, 71);
        String readBases = refBaseStr.substring(20, scStartRefIndex) + scBases;
        int scReadIndex = 31;

        AltRead altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, false);

        assertNotNull(altRead);
        assertEquals("G", altRead.Ref);
        assertEquals("G" + insertBases, altRead.Alt);

        // then on the left
        scStartRefIndex = 0;
        scBases = refBaseStr.substring(scStartRefIndex, 20) + insertBases;
        readBases = scBases + refBaseStr.substring(20, 51);
        scReadIndex = 0;

        altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, true);

        assertNotNull(altRead);
        assertEquals("T", altRead.Ref);
        assertEquals("T" + insertBases, altRead.Alt);
    }

    @Test
    public void testSnvBeforeInsert()
    {
        // create an SNV directly before a 1-base insert - confirm impact on read contexts and how reads are handled

        String refBases = REF_BASES_200.substring(0, 100) + TEST_LEFT_FLANK + "ACGTTCCAACCTTGCA" + REF_BASES_200.substring(0, 100);
        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        //            110               120         130
        // 0123456789 0123456 7      8 9012345  6789012345
        // test flank ACGTTCC A>G insT ACCTTGCA AAAAAGGGGG

        int position = 117;
        SimpleVariant variant = new SimpleVariant(CHR_1, position, "A", "G");

        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 0, 200);

        RefContextCache refContextCache = new RefContextCache(TEST_CONFIG, Collections.emptyList(), Collections.emptyList());

        RefContextConsumer refContextConsumer = new RefContextConsumer(TEST_CONFIG, region, refSequence, refContextCache, Collections.emptyList());

        // send through a read with the SNV, then another with the SNV also immediately followed by an insert
        String readBases = refBases.substring(100, 117) + variant.alt() + refBases.substring(118, 135);

        String cigar = buildCigarString(readBases.length());

        // the read's alignment start with the first base of the read context
        SAMRecord read = buildSamRecord(100, cigar, readBases);
        read.setMappingQuality(60);

        refContextConsumer.processRead(read);
        refContextConsumer.processRead(read); // repeat to add support

        readBases = refBases.substring(100, 117) + variant.alt() + "T" + refBases.substring(118, 135);

        cigar ="18M1I17M";

        read = buildSamRecord(100, cigar, readBases);
        read.setMappingQuality(60);

        refContextConsumer.processRead(read);
        refContextConsumer.processRead(read);

        List<AltContext> altContexts = refContextCache.altContexts();
        assertEquals(2, altContexts.size());

        AltContext snv = altContexts.stream().filter(x -> x.Ref.equals(variant.ref()) && x.Alt.equals(variant.alt())).findFirst().orElse(null);
        assertNotNull(snv);

        AltContext snvInsert = altContexts.stream().filter(x -> x.Ref.equals(variant.ref()) && x.Alt.equals("GT")).findFirst().orElse(null);
        assertNotNull(snvInsert);
    }
}
