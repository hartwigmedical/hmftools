package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SliceRegionTest
{
    private static final SliceConfig CONFIG = new SliceConfig();

    private final ReadCache mReadCache;
    private final SliceWriter mSliceWriter;

    private static final String READ_BASES = "";
    private static final String READ_CIGAR = "100M";

    public SliceRegionTest()
    {
        mReadCache = new ReadCache(CONFIG);
        mSliceWriter = new SliceWriter(CONFIG);
    }

    @Test
    public void testRegionSlicer()
    {
        ChrBaseRegion sliceRegion = new ChrBaseRegion(CHR_1, 1000, 2000);
        RegionBamSlicer regionSlicer = new RegionBamSlicer(sliceRegion, CONFIG, mReadCache, mSliceWriter);

        String readId1 = "READ_01";


        SAMRecord read1 = createSamRecord(readId1, CHR_1, 1000, CHR_1, 1300,false, null);

        regionSlicer.processSamRecord(read1);

        assertEquals(1, regionSlicer.readGroupMap().size());
        assertTrue(mReadCache.chrRemotePositions().isEmpty());

        SAMRecord mate1 = createSamRecord(readId1, CHR_1, 1300, CHR_1, 1000,false, null);

        regionSlicer.processSamRecord(mate1);

        assertTrue(regionSlicer.readGroupMap().isEmpty());
        assertTrue(mReadCache.chrRemotePositions().isEmpty());
    }

    private static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String mateChr, int mateStart,
            boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        return SamRecordTestUtils.createSamRecord(
                readId, chrStr, readStart, READ_BASES, READ_CIGAR, mateChr, mateStart, false,
                false, null, false, READ_CIGAR);
    }

}
