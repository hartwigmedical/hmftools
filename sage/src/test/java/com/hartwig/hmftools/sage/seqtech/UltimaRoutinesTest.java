package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.cloneSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.sage.common.TestUtils.setUltimaSequencing;
import static com.hartwig.hmftools.sage.common.VariantUtils.TEST_LEFT_FLANK;
import static com.hartwig.hmftools.sage.filter.FilterConfig.ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.candidate.ReadContextCandidate;
import com.hartwig.hmftools.sage.candidate.RefContextCache;
import com.hartwig.hmftools.sage.candidate.RefContextConsumer;
import com.hartwig.hmftools.sage.common.RefSequence;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaRoutinesTest
{
    public UltimaRoutinesTest()
    {
        setUltimaSequencing();
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }


    @Test
    public void testLowQualBasesInCore()
    {
        // low-qual bases in core invalidates candidate evidence unless the variant has a long homopolymer
        String refBases = REF_BASES_200.substring(0, 10) + TEST_LEFT_FLANK + "AACTT" + REF_BASES_200.substring(0, 10);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());

        int position = 23;
        SimpleVariant variant = new SimpleVariant(CHR_1, position, "C", "G");

        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);

        RefContextCache refContextCache = new RefContextCache(TEST_CONFIG, Collections.emptyList(), Collections.emptyList());

        RefContextConsumer refContextConsumer = new RefContextConsumer(TEST_CONFIG, region, refSequence, refContextCache, Collections.emptyList());

        // first a valid read without low-qual bases
        String readBases = refBases.substring(0, 22) + variant.alt() + refBases.substring(23);

        String cigar = buildCigarString(readBases.length());

        // the read's alignment start with the first base of the read context
        SAMRecord read = buildSamRecord(1, cigar, readBases);

        refContextConsumer.processRead(read);
        read = cloneSamRecord(read, READ_ID_GENERATOR.nextId());
        refContextConsumer.processRead(read);

        List<ReadContextCandidate> altContexts = refContextCache.altCandidates();
        assertEquals(1, altContexts.size());

        assertTrue(altContexts.stream().anyMatch(x -> x.readContext().variant().matches(variant)));

        // repeat with low-qual bases
        refContextCache.clear();

        String lowCoreQuals = "20=10-29";
        read.setAttribute(ULT_QUAL_TAG, lowCoreQuals);
        refContextConsumer.processRead(read);

        read = cloneSamRecord(read, READ_ID_GENERATOR.nextId());
        refContextConsumer.processRead(read);

        altContexts = refContextCache.altCandidates();
        assertTrue(altContexts.isEmpty());

        // now a high-qual read followed by a low-qual, still passes the min raw alt threshold
        refContextCache.clear();

        read = buildSamRecord(1, cigar, readBases);
        refContextConsumer.processRead(read);

        read = cloneSamRecord(read, READ_ID_GENERATOR.nextId());
        read.setAttribute(ULT_QUAL_TAG, lowCoreQuals);
        refContextConsumer.processRead(read);

        altContexts = refContextCache.altCandidates();
        assertEquals(1, altContexts.size());

        assertTrue(altContexts.stream().anyMatch(x -> x.readContext().variant().matches(variant)));

        // test again for a variant with a long HP repeat
        refContextCache.clear();

        // max repeat will be the GA, but the 9xTs will cause the low-qual check to be skipped
        refBases = REF_BASES_200.substring(0, 10) + TEST_LEFT_FLANK + "AACGAGAGATTTT" + "TTTTTACGTA" + REF_BASES_200.substring(0, 10);
        refSequence = new RefSequence(1, refBases.getBytes());

        refContextConsumer = new RefContextConsumer(TEST_CONFIG, region, refSequence, refContextCache, Collections.emptyList());

        readBases = refBases.substring(0, 22) + variant.alt() + refBases.substring(23);

        cigar = buildCigarString(readBases.length());

        read = buildSamRecord(1, cigar, readBases);

        lowCoreQuals = "10=25-34";
        read.setAttribute(ULT_QUAL_TAG, lowCoreQuals);
        refContextConsumer.processRead(read);

        read = cloneSamRecord(read, READ_ID_GENERATOR.nextId());
        read.setAttribute(ULT_QUAL_TAG, lowCoreQuals);
        refContextConsumer.processRead(read);

        altContexts = refContextCache.altCandidates();
        assertEquals(1, altContexts.size());

        assertTrue(altContexts.stream().anyMatch(x -> x.readContext().variant().matches(variant)));
    }
}
