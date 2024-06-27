package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;

import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class CoverageTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MetricsConfig mConfig;

    protected static final String TEST_READ_BASES = MockRefGenome.generateRandomBases(10);
    protected static final String TEST_CIGAR = "10M";

    public CoverageTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mConfig = new MetricsConfig(10);
    }

    @Test
    public void testReadFilterTypes()
    {
        BaseCoverage baseCoverage = new BaseCoverage(mConfig, 1, 200, Collections.emptyList());

        // add 1 of each filter type:

        // unfiltered
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);

        baseCoverage.processRead(read, null, false);

        int readLength = read.getReadLength();

        // low map qual
        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);
        read.setMappingQuality(mConfig.MapQualityThreshold - 1);

        baseCoverage.processRead(read, null, false);

        // unmapped
        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, "*", 0,
                false, false, null);
        read.setMateUnmappedFlag(true);

        baseCoverage.processRead(read, null, false);

        // duplicate x 2 and a consensus read
        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);
        read.setDuplicateReadFlag(true);

        baseCoverage.processRead(read, null, false);
        baseCoverage.processRead(read, null, false);

        // consensus
        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);
        addConsensusReadAttribute(read, 2, 1, UmiReadType.SINGLE);

        baseCoverage.processRead(read, null, true);

        // low base qual
        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);

        for(int i = 0; i < read.getBaseQualities().length; ++i)
        {
            read.getBaseQualities()[i] = (byte) (mConfig.BaseQualityThreshold - 1);
        }

        baseCoverage.processRead(read, null, false);

        CoverageMetrics metrics = baseCoverage.createMetrics();

        assertEquals(readLength * 2, metrics.FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(readLength, metrics.FilterTypeCounts[FilterType.LOW_MAP_QUAL.ordinal()]);
        assertEquals(readLength, metrics.FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
        assertEquals(readLength, metrics.FilterTypeCounts[FilterType.LOW_BASE_QUAL.ordinal()]);

        baseCoverage.clear();

        // exceeding max coverage - bases 30-39 will now have 30x coverage
        String testBases = MockRefGenome.generateRandomBases(30);
        String testCigar = "5S20M5S";

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 10, testBases, testCigar, CHR_1, 100,
                false, false, null);

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 15, testBases, testCigar, CHR_1, 100,
                false, false, null);

        SAMRecord read3 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, testBases, testCigar, CHR_1, 100,
                false, false, null);

        for(int i = 0; i < 10; ++i)
        {
            baseCoverage.processRead(read1, null, false);
            baseCoverage.processRead(read2, null, false);
            baseCoverage.processRead(read3, null, false);
        }

        metrics = baseCoverage.createMetrics();

        for(int i = 0; i < baseCoverage.baseDepth().length; ++i)
        {
            assertTrue(baseCoverage.baseDepth()[i] <= mConfig.MaxCoverage);
        }

        assertEquals(300, metrics.FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(300, metrics.FilterTypeCounts[FilterType.MAX_COVERAGE.ordinal()]);
    }

    @Test
    public void testOverlappingReadCoverage()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        BamReader bamReader = new BamReader(
                new ChrBaseRegion(CHR_1, 1, 500), mConfig, null, null, combinedStats);

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 25,
                false, false, null);

        bamReader.processRead(read1);

        // now its mate
        SAMRecord mate1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 25, TEST_READ_BASES, TEST_CIGAR, CHR_1, 20,
                true, false, null);

        bamReader.processRead(mate1);

        BaseCoverage baseCoverage = bamReader.baseCoverage();
        CoverageMetrics metrics = baseCoverage.createMetrics();

        assertEquals(15, metrics.FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);

        assertTrue(bamReader.readGroupMap().isEmpty());

        for(int i = 0; i < baseCoverage.baseDepth().length; ++i)
        {
            if(i >= 19 && i <= 33)
                assertTrue(baseCoverage.baseDepth()[i] == 1);
            else
                assertTrue(baseCoverage.baseDepth()[i] == 0);
        }

        // a more complicate example for overlapping alignments..
        baseCoverage.clear();

        String testBases = MockRefGenome.generateRandomBases(100);

        // actual bases don't need to match

        String readCigar1 = "5S10M10D10M5I10M5S";
        String mateCigar1 = "5S10M5I10M10N10M5S";
        String suppCigar1 = "5S10M10N10M5S";

        read1 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 10, testBases.substring(0, 45), readCigar1, CHR_1, 30,
                false, false, new SupplementaryReadData(CHR_1, 40, SUPP_POS_STRAND, suppCigar1, 60));

        bamReader.processRead(read1);

        // now its mate
        mate1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 30, testBases.substring(0, 45), mateCigar1, CHR_1, 10,
                true, false, null);

        bamReader.processRead(mate1);

        SAMRecord supp1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 40, testBases.substring(0, 45), suppCigar1, CHR_1, 30,
                false, true, new SupplementaryReadData(CHR_1, 10, SUPP_POS_STRAND, readCigar1, 60));

        bamReader.processRead(supp1);

        for(int i = 0; i < baseCoverage.baseDepth().length; ++i)
        {
            if((i >= 9 && i <= 18) || (i >= 29 && i <= 48) || (i >= 59 && i <= 68))
                assertTrue(baseCoverage.baseDepth()[i] == 1);
            else
                assertTrue(baseCoverage.baseDepth()[i] == 0);
        }

        assertTrue(bamReader.readGroupMap().isEmpty());

        baseCoverage.clear();

        // supplementary overlaps but not the mate - still checks overlapping bases
        String testCigar = "20M";
        read1 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 10, testBases.substring(0, 20), testCigar, CHR_1, 30,
                false, false, new SupplementaryReadData(CHR_1, 20, SUPP_POS_STRAND, testCigar, 60));

        bamReader.processRead(read1);

        // now its mate
        mate1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 50, testBases.substring(50, 70), testCigar, CHR_1, 10,
                true, false, null);

        bamReader.processRead(mate1);

        supp1 = SamRecordTestUtils.createSamRecord(
                read1.getReadName(), CHR_1, 20, testBases.substring(10, 30), testCigar, CHR_1, 10,
                false, true, new SupplementaryReadData(CHR_1, 10, SUPP_POS_STRAND, testCigar, 60));

        bamReader.processRead(supp1);

        for(int i = 0; i < baseCoverage.baseDepth().length; ++i)
        {
            if((i >= 9 && i <= 38) || (i >= 49 && i <= 68))
                assertTrue(baseCoverage.baseDepth()[i] == 1);
            else
                assertTrue(baseCoverage.baseDepth()[i] == 0);
        }

        assertTrue(bamReader.readGroupMap().isEmpty());
    }
}
