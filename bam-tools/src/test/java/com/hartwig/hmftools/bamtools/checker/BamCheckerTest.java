package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.bam.BamReadLite;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class BamCheckerTest
{
    private static final ReadIdGenerator READ_ID_GENERATOR = new ReadIdGenerator();
    private static final String TEST_READ_BASES = "ACGTACGTACGTACGTACGT";
    private static final String TEST_CIGAR = "20M";

    @Test
    public void testBasicFragment()
    {
        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, null);

        Fragment fragment = new Fragment(read1);

        assertTrue(!fragment.hasPrimaryInfo());
        assertTrue(!fragment.isComplete());
        assertTrue(fragment.extractCompleteReads().isEmpty());

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 300, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, false, null);
        SamRecordTestUtils.flipFirstInPair(read2);

        fragment.addRead(read2);

        assertTrue(fragment.hasPrimaryInfo());
        assertTrue(fragment.isComplete());

        List<SAMRecord> completeReads = fragment.extractCompleteReads();
        assertEquals(2, completeReads.size());
        assertTrue(fragment.reads().isEmpty());
    }

    @Test
    public void testFragmentWithSupplementaries()
    {
        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, new SupplementaryReadData(CHR_1, 200, SUPP_POS_STRAND, TEST_CIGAR, 60));

        SAMRecord supp1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 200, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, true, new SupplementaryReadData(CHR_1, 100, SUPP_POS_STRAND, TEST_CIGAR, 60));

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 300, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, false, new SupplementaryReadData(CHR_1, 400, SUPP_POS_STRAND, TEST_CIGAR, 60));
        SamRecordTestUtils.flipFirstInPair(read2);

        SAMRecord supp2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 400, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, true, new SupplementaryReadData(CHR_1, 300, SUPP_POS_STRAND, TEST_CIGAR, 60));
        SamRecordTestUtils.flipFirstInPair(supp2);

        // test each distinct combination of reads being received:

        // test 1: R1, both supps then R2
        Fragment fragment = new Fragment(read1);

        assertTrue(!fragment.hasPrimaryInfo());
        assertTrue(!fragment.isComplete());

        assertEquals(1, fragment.expectedSupplementaryCount());
        assertEquals(0, fragment.receivedSupplementaryCount());

        fragment.addRead(supp1);
        fragment.addRead(supp2);

        assertEquals(1, fragment.expectedSupplementaryCount());
        assertEquals(2, fragment.receivedSupplementaryCount());

        assertTrue(!fragment.hasPrimaryInfo());
        assertTrue(!fragment.isComplete());

        fragment.addRead(read2);

        assertTrue(fragment.hasPrimaryInfo());
        assertTrue(fragment.isComplete());

        List<SAMRecord> completeReads = fragment.extractCompleteReads();
        assertEquals(4, completeReads.size());
        assertTrue(fragment.reads().isEmpty());

        // test 2: both supps then both primaries
        fragment = new Fragment(supp1);
        fragment.addRead(supp2);

        assertTrue(!fragment.hasPrimaryInfo());
        assertTrue(!fragment.isComplete());

        assertEquals(0, fragment.expectedSupplementaryCount());
        assertEquals(2, fragment.receivedSupplementaryCount());

        fragment.addRead(read2);

        assertTrue(!fragment.hasPrimaryInfo());
        assertTrue(!fragment.isComplete());

        assertEquals(1, fragment.expectedSupplementaryCount());

        fragment.addRead(read1);

        assertTrue(fragment.hasPrimaryInfo());
        assertTrue(fragment.isComplete());

        completeReads = fragment.extractCompleteReads();
        assertEquals(4, completeReads.size());
        assertTrue(fragment.reads().isEmpty());

        // test 3: supps coming last
        fragment = new Fragment(read2);
        fragment.addRead(supp1);

        assertTrue(!fragment.hasPrimaryInfo());
        assertTrue(!fragment.isComplete());

        assertEquals(1, fragment.expectedSupplementaryCount());
        assertEquals(1, fragment.receivedSupplementaryCount());

        fragment.addRead(read1);

        assertTrue(fragment.hasPrimaryInfo());
        assertTrue(!fragment.isComplete());

        assertEquals(2, fragment.expectedSupplementaryCount());

        completeReads = fragment.extractCompleteReads();
        assertEquals(3, completeReads.size());
        assertTrue(fragment.reads().isEmpty());

        fragment.addRead(supp2);

        assertTrue(fragment.isComplete());

        completeReads = fragment.extractCompleteReads();
        assertEquals(1, completeReads.size());
        assertTrue(fragment.reads().isEmpty());
    }

    @Test
    public void testFragmentCache()
    {
        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, null);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 300, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, false, null);
        SamRecordTestUtils.flipFirstInPair(read2);


        FragmentCache fragmentCache = new FragmentCache(null);

        // test 1: R1 and R2, no supplementaries
        Fragment fragment = new Fragment(read1);

        List<SAMRecord> completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));
        assertTrue(completeReads.isEmpty());
        assertEquals(1, fragmentCache.fragmentCount());

        fragment = new Fragment(read2);

        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));
        assertEquals(2, completeReads.size());
        assertEquals(read1, completeReads.get(0));
        assertEquals(read2, completeReads.get(1));
        assertEquals(0, fragmentCache.fragmentCount());

        read1 = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, new SupplementaryReadData(CHR_1, 200, SUPP_POS_STRAND, TEST_CIGAR, 60));

        SAMRecord supp1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 200, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, true, new SupplementaryReadData(CHR_1, 100, SUPP_POS_STRAND, TEST_CIGAR, 60));

        read2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 300, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, false, new SupplementaryReadData(CHR_1, 400, SUPP_POS_STRAND, TEST_CIGAR, 60));
        SamRecordTestUtils.flipFirstInPair(read2);

        SAMRecord supp2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 400, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, true, new SupplementaryReadData(CHR_1, 300, SUPP_POS_STRAND, TEST_CIGAR, 60));
        SamRecordTestUtils.flipFirstInPair(supp2);

        // test 2: R1, both supps then R2
        fragment = new Fragment(read1);

        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));
        assertTrue(completeReads.isEmpty());
        assertEquals(1, fragmentCache.fragmentCount());

        fragment = new Fragment(supp2);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertTrue(completeReads.isEmpty());
        assertEquals(1, fragmentCache.fragmentCount());
        assertEquals(2, fragmentCache.readCount());

        fragment = new Fragment(supp1);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertTrue(completeReads.isEmpty());
        assertEquals(1, fragmentCache.fragmentCount());
        assertEquals(3, fragmentCache.readCount());

        fragment = new Fragment(read2);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertEquals(4, completeReads.size());
        assertEquals(0, fragmentCache.fragmentCount());
        assertEquals(0, fragmentCache.readCount());

        // test 3: primaries then supplementaries
        fragment = new Fragment(read2);

        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));
        assertTrue(completeReads.isEmpty());
        assertEquals(1, fragmentCache.fragmentCount());

        fragment = new Fragment(read1);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertEquals(2, completeReads.size());
        assertEquals(1, fragmentCache.fragmentCount());
        assertEquals(0, fragmentCache.readCount());

        fragment = new Fragment(supp1);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertEquals(1, completeReads.size());
        assertEquals(supp1, completeReads.get(0));
        assertEquals(1, fragmentCache.fragmentCount());
        assertEquals(0, fragmentCache.readCount());

        fragment = new Fragment(supp2);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertEquals(1, completeReads.size());
        assertEquals(supp2, completeReads.get(0));
        assertEquals(0, fragmentCache.fragmentCount());
        assertEquals(0, fragmentCache.readCount());

        // test 4: supp is cached and then pass through fragment which had complete R1 and R2, just to set mate info
        fragment = new Fragment(supp1);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertTrue(completeReads.isEmpty());
        assertEquals(1, fragmentCache.fragmentCount());
        assertEquals(1, fragmentCache.readCount());

        fragment = new Fragment(supp2);
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));

        assertTrue(completeReads.isEmpty());
        assertEquals(1, fragmentCache.fragmentCount());
        assertEquals(2, fragmentCache.readCount());

        fragment = new Fragment(read1);
        fragment.addRead(read2);
        assertTrue(fragment.hasPrimaryInfo());
        assertEquals(2, fragment.expectedSupplementaryCount());

        completeReads = fragment.extractCompleteReads();
        assertEquals(2, completeReads.size());
        assertEquals(0, fragment.readCount());

        // now the supps will have their mate CIGARs set and be processed
        completeReads = fragmentCache.handleIncompleteFragments(List.of(fragment));
        assertEquals(2, completeReads.size());
        assertEquals(0, fragmentCache.fragmentCount());
        assertEquals(0, fragmentCache.readCount());
    }

    @Test
    public void testReadLite()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, new SupplementaryReadData(CHR_1, 200, SUPP_POS_STRAND, TEST_CIGAR, 60));

        read.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_CIGAR);

        BamReadLite readLite = new BamReadLite(read, true);

        SAMRecord restoredRead = BamReadLite.from(readLite, null, read.getReadName());
        assertEquals(read.getReadName(), restoredRead.getReadName());
        assertEquals(read.getAlignmentStart(), restoredRead.getAlignmentStart());
        assertEquals(read.getReferenceName(), restoredRead.getReferenceName());
        assertEquals(read.getCigarString(), restoredRead.getCigarString());
        assertEquals(read.getInferredInsertSize(), restoredRead.getInferredInsertSize());
        assertEquals(read.getMateReferenceName(), restoredRead.getMateReferenceName());
        assertEquals(read.getMateAlignmentStart(), restoredRead.getMateAlignmentStart());

        for(SAMRecord.SAMTagAndValue attribute : read.getAttributes())
        {
            Object otherAttribute = restoredRead.getAttribute(attribute.tag);
            assertNotNull(otherAttribute);
            assertEquals(otherAttribute, read.getAttribute(attribute.tag));
        }
    }
}
