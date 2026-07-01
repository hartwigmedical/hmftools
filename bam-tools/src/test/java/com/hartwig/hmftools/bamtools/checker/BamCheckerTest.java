package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.bamtools.checker.Fragment.convertHardClips;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.READ_GROUP_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.S;

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

    public BamCheckerTest()
    {
        CheckConfig.Params.MinAlignmentScore = 10;
    }

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
    public void testUnmappedReads()
    {
        // first R1 requiring unmapping
        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, new SupplementaryReadData(CHR_1, 200, SUPP_POS_STRAND, TEST_CIGAR, 60));

        int lowAlignmentScore = CheckConfig.Params.MinAlignmentScore - 1;
        read1.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, lowAlignmentScore);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 300, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, false, new SupplementaryReadData(CHR_1, 400, SUPP_POS_STRAND, TEST_CIGAR, 60));
        SamRecordTestUtils.flipFirstInPair(read2);

        SAMRecord supp1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 200, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, true, new SupplementaryReadData(CHR_1, 100, SUPP_POS_STRAND, TEST_CIGAR, 60));

        SAMRecord supp2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 400, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, true, new SupplementaryReadData(CHR_1, 300, SUPP_POS_STRAND, TEST_CIGAR, 60));
        SamRecordTestUtils.flipFirstInPair(supp2);

        Fragment fragment = new Fragment(read1);
        fragment.addRead(read2);
        fragment.addRead(supp1);
        fragment.addRead(supp2);
        List<SAMRecord> reads = fragment.extractCompleteReads();

        assertEquals(1, fragment.expectedSupplementaryCount());
        assertEquals(1, fragment.receivedSupplementaryCount());
        assertEquals(3, reads.size());

        assertFalse(reads.contains(supp1));

        assertTrue(read1.getReadUnmappedFlag());
        assertEquals(NO_CIGAR, read1.getCigarString());
        assertEquals(0, read1.getInferredInsertSize());
        assertNull(read1.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertEquals(read2.getAlignmentStart(), read1.getAlignmentStart());

        assertTrue(read2.getMateUnmappedFlag());

        // test again with the supp arriving first
        fragment = new Fragment(supp1);
        fragment.addRead(read2);
        fragment.addRead(read1);
        fragment.addRead(supp2);
        reads = fragment.extractCompleteReads();

        assertEquals(3, reads.size());

        assertFalse(reads.contains(supp1));
        assertTrue(read1.getReadUnmappedFlag());
        assertNull(read1.getAttribute(SUPPLEMENTARY_ATTRIBUTE));

        // for the second read
        read1 = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, new SupplementaryReadData(CHR_1, 200, SUPP_POS_STRAND, TEST_CIGAR, 60));

        read2.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, lowAlignmentScore);

        fragment = new Fragment(supp2);
        fragment.addRead(read2);
        fragment.addRead(read1);
        fragment.addRead(supp1);
        reads = fragment.extractCompleteReads();

        assertEquals(3, reads.size());
        assertEquals(1, fragment.expectedSupplementaryCount());
        assertEquals(1, fragment.receivedSupplementaryCount());

        assertTrue(reads.contains(supp1));
        assertFalse(reads.contains(supp2));
        assertTrue(read2.getReadUnmappedFlag());
        assertNull(read2.getAttribute(SUPPLEMENTARY_ATTRIBUTE));

        // for both reads
        read1.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, lowAlignmentScore);

        read2 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 300, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, false, new SupplementaryReadData(CHR_1, 400, SUPP_POS_STRAND, TEST_CIGAR, 60));
        SamRecordTestUtils.flipFirstInPair(read2);
        read2.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, lowAlignmentScore);

        fragment = new Fragment(read1);
        fragment.addRead(read2);
        fragment.addRead(supp1);
        fragment.addRead(supp2);
        reads = fragment.extractCompleteReads();

        assertEquals(2, reads.size());
        assertEquals(0, fragment.expectedSupplementaryCount());
        assertEquals(0, fragment.receivedSupplementaryCount());

        assertFalse(reads.contains(supp1));
        assertFalse(reads.contains(supp2));
        assertTrue(read1.getReadUnmappedFlag());
        assertTrue(read2.getReadUnmappedFlag());
        assertNull(read1.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertNull(read2.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
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
        // assertTrue(completeReads.contains(read1));
        //assertTrue(completeReads.contains(read2));
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

        String readGroupId = "RG_ID_001";

        read.setAttribute(READ_GROUP_ATTRIBUTE, readGroupId);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, 0);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 20);
        read.setAttribute(XS_ATTRIBUTE, 30);

        BamReadLite readLite = new BamReadLite(read, true);

        SAMRecord restoredRead = BamReadLite.from(readLite, null, read.getReadName(), readGroupId);
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

        // test again with no additional attributes

        read = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 300,
                false, false, null);

        read.setAttribute(READ_GROUP_ATTRIBUTE, readGroupId);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, 20);

        readLite = new BamReadLite(read, true);

        restoredRead = BamReadLite.from(readLite, null, read.getReadName(), readGroupId);
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

    @Test
    public void testHardClippedSupplementaries()
    {
        String suppReadBases = "A".repeat(30);
        String primaryReadBases = "A".repeat(60);
        byte[] primaryReadQuals = buildBaseQuals(primaryReadBases.length(), 30);

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, suppReadBases, "10H30M20H", CHR_1, 300,
                false, false, new SupplementaryReadData(CHR_1, 200, SUPP_POS_STRAND, TEST_CIGAR, 60));

        convertHardClips(read, primaryReadBases.getBytes(), primaryReadQuals, false);

        assertEquals(3, read.getCigar().getCigarElements().size());
        assertTrue(read.getCigar().getCigarElements().get(0).getOperator() == S);
        assertTrue(read.getCigar().getCigarElements().get(2).getOperator() == S);

        assertEquals(primaryReadBases.getBytes().length, read.getReadBases().length);
    }
}
