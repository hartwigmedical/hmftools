package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Queue;

import com.google.common.collect.Queues;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
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
        mSliceWriter = new SliceWriter(CONFIG);
        mReadCache = new ReadCache(mSliceWriter);
    }

    /*
    // TODO: fix this
    @Test
    public void testRegionSlicer()
    {
        ChrBaseRegion sliceRegion = new ChrBaseRegion(CHR_1, 1000, 2000);
        ChrBaseRegion sliceRegion2 = new ChrBaseRegion(CHR_1, 2050, 2500);
        ChrBaseRegion sliceRegion3 = new ChrBaseRegion(CHR_1, 2550, 3000);

        CONFIG.SliceRegions.addRegion(sliceRegion);
        CONFIG.SliceRegions.addRegion(sliceRegion2);
        CONFIG.SliceRegions.addRegion(sliceRegion3);

        Queue<ChrBaseRegion> regions = Queues.newArrayDeque();
        RegionBamSlicer regionSlicer = new RegionBamSlicer(regions, CONFIG, mReadCache, mSliceWriter);
        regionSlicer.setCurrentRegion(sliceRegion);

        ReadIdGenerator readIdGenerator = new ReadIdGenerator();

        // before the region
        SAMRecord read = createSamRecord(readIdGenerator.nextId(), CHR_1, 850, CHR_1, 1300,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(0, regionSlicer.readsProcessed());

        read = createSamRecord(readIdGenerator.nextId(), CHR_1, 950, CHR_1, 1300,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(1, regionSlicer.readsProcessed());

        // mate is within the region
        assertTrue(mReadCache.chrRemotePositions().isEmpty());

        read = createSamRecord(readIdGenerator.nextId(), CHR_1, 1300, CHR_1, 1000,false, null);
        regionSlicer.processSamRecord(read);

        assertEquals(2, regionSlicer.readsProcessed());
        assertTrue(mReadCache.chrRemotePositions().isEmpty());

        // mate just outside region so not captured by initial slice
        read = createSamRecord(readIdGenerator.nextId(), CHR_1, 1200, CHR_1, 800,false, null);
        regionSlicer.processSamRecord(read);

        assertEquals(3, regionSlicer.readsProcessed());
        assertEquals(1, mReadCache.chrRemotePositions().size());

        // test ignoring reads overlapped by an earlier region
        regionSlicer.setCurrentRegion(sliceRegion2);

        read = createSamRecord(readIdGenerator.nextId(), CHR_1, 1980, CHR_2, 1000,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(3, regionSlicer.readsProcessed());
        assertEquals(1, mReadCache.chrRemotePositions().size());

        read = createSamRecord(readIdGenerator.nextId(), CHR_1, 2020, CHR_2, 1000,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(4, regionSlicer.readsProcessed());
        assertEquals(2, mReadCache.chrRemotePositions().size());

        regionSlicer.setCurrentRegion(sliceRegion3);

        read = createSamRecord(readIdGenerator.nextId(), CHR_1, 2480, CHR_2, 1000,false, null);
        regionSlicer.processSamRecord(read);
        assertEquals(4, regionSlicer.readsProcessed());
        assertEquals(2, mReadCache.chrRemotePositions().size());

        // mate in later region, supplementary on diff chromosome
        read = createSamRecord(
                readIdGenerator.nextId(), CHR_1, 2600, CHR_1, 3500,
                false, new SupplementaryReadData(CHR_2, 10000, SUPP_NEG_STRAND, READ_CIGAR, 60));
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        regionSlicer.processSamRecord(read);

        assertEquals(2, mReadCache.chrRemotePositions().size());
        assertEquals(2, mReadCache.chrRemotePositions().get(CHR_1).size());
        assertEquals(2, mReadCache.chrRemotePositions().get(CHR_2).size());

    }*/

    private static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String mateChr, int mateStart,
            boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        return SamRecordTestUtils.createSamRecord(
                readId, chrStr, readStart, READ_BASES, READ_CIGAR, mateChr, mateStart, false,
                isSupplementary, suppAlignment, false, READ_CIGAR);
    }

}
