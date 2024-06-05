package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.Queue;

import com.google.common.collect.Queues;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
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

        Queue<ChrBaseRegion> regions = Queues.newArrayDeque();
        RegionBamSlicer regionSlicer = new RegionBamSlicer(regions, CONFIG, mReadCache, mSliceWriter);
        regionSlicer.setCurrentRegion(sliceRegion);

        String readId = "READ_01";

        SAMRecord read = createSamRecord(readId, CHR_1, 1000, CHR_1, 1300,false, null);

        regionSlicer.processSamRecord(read);

        assertTrue(mReadCache.chrRemotePositions().isEmpty());

        read = createSamRecord(readId, CHR_1, 1300, CHR_1, 1000,false, null);

        regionSlicer.processSamRecord(read);

        assertTrue(mReadCache.chrRemotePositions().isEmpty());

        // mate just outside region so not captured by initial slice
        read = createSamRecord(readId, CHR_1, 1200, CHR_1, 990,false, null);

        regionSlicer.processSamRecord(read);

        assertEquals(1, mReadCache.chrRemotePositions().size());

        mReadCache.chrRemotePositions().clear();

        // mate in later region, supplementary on diff chromosome
        read = createSamRecord(
                readId, CHR_1, 1100, CHR_1, 2010,
                false, new SupplementaryReadData(CHR_2, 1020, SUPP_NEG_STRAND, READ_CIGAR, 60));

        regionSlicer.processSamRecord(read);

        assertEquals(2, mReadCache.chrRemotePositions().size());

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
