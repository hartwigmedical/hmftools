package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.metrics.CoverageTest.TEST_CIGAR;
import static com.hartwig.hmftools.bamtools.metrics.CoverageTest.TEST_READ_BASES;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MetricsTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MetricsConfig mConfig;

    public MetricsTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mConfig = new MetricsConfig(10);
    }

    @Test
    public void testReadCounts()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        BamReader bamReader = new BamReader(new ChrBaseRegion(
                CHR_1, 1, 1000), mConfig, null, null, combinedStats);

        // a primary, non-consensus
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);

        bamReader.processRead(read);

        // 2x duplicates and consensus
        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);

        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);

        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, false, null);

        addConsensusReadAttribute(read, 2, 1, UmiReadType.SINGLE);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(3, combinedStats.readCounts().Total);
        assertEquals(1, combinedStats.readCounts().Duplicates);
        assertEquals(3, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
    }

    @Test
    public void testDualStrandSupplementaries()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        BamReader bamReader = new BamReader(new ChrBaseRegion(
                CHR_1, 1, 1000), mConfig, null, null, combinedStats);

        // supplementary with a duplicate
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, true, null);

        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                true, true, null);

        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, true, null);

        addConsensusReadAttribute(read, 2, 1, UmiReadType.DUAL);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(2, combinedStats.readCounts().Total);
        assertEquals(1, combinedStats.readCounts().Duplicates);
        assertEquals(1, combinedStats.readCounts().DualStrand);
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(2, combinedStats.flagStats().passCount(FlagStatType.SUPPLEMENTARY));
    }
}
