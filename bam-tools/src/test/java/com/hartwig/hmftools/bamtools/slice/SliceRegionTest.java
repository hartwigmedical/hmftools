package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.ALIGNMENTS_DELIM;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_4;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SliceRegionTest
{
    private static final SliceConfig CONFIG = new SliceConfig();
    private static final ReadIdGenerator READ_ID_GENERATOR = new ReadIdGenerator();

    private final ReadCache mReadCache;
    private final SliceWriter mSliceWriter;

    private static final String READ_BASES = "";
    private static final String READ_CIGAR = "100M";

    public SliceRegionTest()
    {
        mSliceWriter = new SliceWriter(CONFIG);
        mReadCache = new ReadCache(mSliceWriter);
    }

    @Test
    public void testRegionSlicerSingleRegion()
    {
        ChrBaseRegion sliceRegion = new ChrBaseRegion(CHR_1, 1000, 2000);

        RegionBamSlicer regionSlicer = new RegionBamSlicer(sliceRegion, CONFIG, mReadCache, null);

        // reads outside the region are ignored
        SAMRecord read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 850, CHR_1, 1300,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(0, regionSlicer.readsProcessed());
        assertEquals(0, mReadCache.fragmentMap().size());

        read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 2001, CHR_1, 1300,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(0, regionSlicer.readsProcessed());
        assertEquals(0, mReadCache.fragmentMap().size());

        // test 1: first read has its mate within the same region
        read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 950, CHR_1, 1300,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(1, mReadCache.fragmentMap().size());

        read = createSamRecord(read.getReadName(), CHR_1, 1300, CHR_1, 950,false, null);
        SamRecordTestUtils.flipFirstInPair(read);
        regionSlicer.processSamRecord(read);
        assertEquals(0, mReadCache.fragmentMap().size());
        assertEquals(2, mSliceWriter.writeCount());

        // test 2: both primaries and their supplementaries in the same region
        SupplementaryReadData suppData = new SupplementaryReadData(CHR_1, 1100, SUPP_POS_STRAND, READ_CIGAR, 60);
        read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 950, CHR_1, 1300,false, suppData);
        regionSlicer.processSamRecord(read);
        assertEquals(1, mReadCache.fragmentMap().size());

        Fragment fragment = mReadCache.fragmentMap().get(read.getReadName());
        assertNotNull(fragment);
        assertEquals(2, fragment.pendingReads().size());

        suppData = new SupplementaryReadData(CHR_1, 950, SUPP_POS_STRAND, READ_CIGAR, 60);
        read = createSamRecord(read.getReadName(), CHR_1, 1100, CHR_1, 1300,true, suppData);
        regionSlicer.processSamRecord(read);
        assertEquals(1, mReadCache.fragmentMap().size());

        SupplementaryReadData suppData1 = new SupplementaryReadData(CHR_1, 1500, SUPP_POS_STRAND, READ_CIGAR, 60);
        SupplementaryReadData suppData2 = new SupplementaryReadData(CHR_1, 1800, SUPP_POS_STRAND, READ_CIGAR, 60);

        // the mate with 2 supps
        read = createSamRecord(read.getReadName(), CHR_1, 1300, CHR_1, 950,false, null);
        SamRecordTestUtils.flipFirstInPair(read);

        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, format("%s%s%s", suppData1.asSamTag(), ALIGNMENTS_DELIM, suppData2.asSamTag()));
        regionSlicer.processSamRecord(read);

        assertEquals(1, mReadCache.fragmentMap().size());
        assertEquals(5, mSliceWriter.writeCount());

        assertEquals(2, fragment.pendingReads().size());

        // supps from second primary - note that the supps also have each other's alignments
        suppData = new SupplementaryReadData(CHR_1, 1300, SUPP_POS_STRAND, READ_CIGAR, 60);
        read = createSamRecord(read.getReadName(), CHR_1, 1500, CHR_1, 950,true, null);
        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, format("%s%s%s", suppData.asSamTag(), ALIGNMENTS_DELIM, suppData2.asSamTag()));
        SamRecordTestUtils.flipFirstInPair(read);

        regionSlicer.processSamRecord(read);
        assertEquals(1, mReadCache.fragmentMap().size());

        read = createSamRecord(read.getReadName(), CHR_1, 1800, CHR_1, 950,true, null);
        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, format("%s%s%s", suppData.asSamTag(), ALIGNMENTS_DELIM, suppData1.asSamTag()));
        SamRecordTestUtils.flipFirstInPair(read);

        regionSlicer.processSamRecord(read);
        assertEquals(0, mReadCache.fragmentMap().size());
    }

    @Test
    public void testRegionSlicerMultiRegion()
    {
        ChrBaseRegion sliceRegion = new ChrBaseRegion(CHR_1, 1000, 2000);
        RegionBamSlicer regionSlicer = new RegionBamSlicer(sliceRegion, CONFIG, mReadCache, null);

        // test 1: all reads are in different regions: primary first chr1, primary second chr 2, supp second chr 3, supp first chr 4
        SupplementaryReadData suppData = new SupplementaryReadData(CHR_4, 1000, SUPP_POS_STRAND, READ_CIGAR, 60);
        SAMRecord read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 1000, CHR_2, 1000, false, suppData);
        regionSlicer.processSamRecord(read);
        assertEquals(1, mReadCache.fragmentMap().size());

        Fragment fragment = mReadCache.fragmentMap().get(read.getReadName());
        assertNotNull(fragment);
        assertEquals(2, fragment.pendingReads().size());

        // the mate
        sliceRegion = new ChrBaseRegion(CHR_2, 1000, 2000);
        regionSlicer = new RegionBamSlicer(sliceRegion, CONFIG, mReadCache, null);

        suppData = new SupplementaryReadData(CHR_3, 1000, SUPP_POS_STRAND, READ_CIGAR, 60);
        read = createSamRecord(read.getReadName(), CHR_2, 1000, CHR_1, 1000, false, suppData);
        SamRecordTestUtils.flipFirstInPair(read);
        regionSlicer.processSamRecord(read);
        assertEquals(1, mReadCache.fragmentMap().size());
        assertEquals(2, mSliceWriter.writeCount());

        assertEquals(2, fragment.pendingReads().size());

        // the mate's supp
        sliceRegion = new ChrBaseRegion(CHR_3, 1000, 2000);
        regionSlicer = new RegionBamSlicer(sliceRegion, CONFIG, mReadCache, null);

        suppData = new SupplementaryReadData(CHR_2, 1000, SUPP_POS_STRAND, READ_CIGAR, 60);
        read = createSamRecord(read.getReadName(), CHR_3, 1000, CHR_1, 1000, true, suppData);
        SamRecordTestUtils.flipFirstInPair(read);
        regionSlicer.processSamRecord(read);
        assertEquals(1, mReadCache.fragmentMap().size());

        assertEquals(1, fragment.pendingReads().size());

        // the first primary's supp
        sliceRegion = new ChrBaseRegion(CHR_4, 1000, 2000);
        regionSlicer = new RegionBamSlicer(sliceRegion, CONFIG, mReadCache, null);

        suppData = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, READ_CIGAR, 60);
        read = createSamRecord(read.getReadName(), CHR_4, 1000, CHR_2, 1000, true, suppData);
        regionSlicer.processSamRecord(read);
        assertEquals(0, mReadCache.fragmentMap().size());

        assertEquals(0, fragment.pendingReads().size());
    }

    @Test
    public void testRegionSlicerUnpaired()
    {
        ChrBaseRegion sliceRegion = new ChrBaseRegion(CHR_1, 1000, 2000);

        RegionBamSlicer regionSlicer = new RegionBamSlicer(sliceRegion, CONFIG, mReadCache, null);

        SupplementaryReadData suppData = new SupplementaryReadData(CHR_1, 1500, SUPP_POS_STRAND, READ_CIGAR, 60);
        SAMRecord read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 1000, NO_CHROMOSOME_NAME, NO_POSITION, false, suppData);
        read.setReadPairedFlag(false);
        regionSlicer.processSamRecord(read);
        assertEquals(1, regionSlicer.readsProcessed());
        assertEquals(1, mReadCache.fragmentMap().size());

        Fragment fragment = mReadCache.fragmentMap().get(read.getReadName());
        assertNotNull(fragment);
        assertTrue(fragment.hasPrimaries());
        assertEquals(1, fragment.pendingReads().size());

        // the supplementary
        suppData = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, READ_CIGAR, 60);
        read = createSamRecord(read.getReadName(), CHR_1, 1500, NO_CHROMOSOME_NAME, NO_POSITION, true, suppData);
        read.setReadPairedFlag(false);
        regionSlicer.processSamRecord(read);
        assertEquals(2, regionSlicer.readsProcessed());
        assertEquals(0, mReadCache.fragmentMap().size());

        assertEquals(0, fragment.pendingReads().size());
    }

    private static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String mateChr, int mateStart,
            boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        return SamRecordTestUtils.createSamRecord(
                readId, chrStr, readStart, READ_BASES, READ_CIGAR, mateChr, mateStart, false,
                isSupplementary, suppAlignment, false, READ_CIGAR);
    }
}
