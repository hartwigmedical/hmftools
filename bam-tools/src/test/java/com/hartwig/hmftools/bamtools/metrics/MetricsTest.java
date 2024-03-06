package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.metrics.CoverageTest.TEST_CIGAR;
import static com.hartwig.hmftools.bamtools.metrics.CoverageTest.TEST_READ_BASES;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.UmiReadType;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.jetbrains.annotations.NotNull;
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
                CHR_1, 1, 1000), mConfig, null, null, combinedStats, null);

        SAMRecord read = getTestRead(false);

        bamReader.processRead(read);

        // 2x duplicates and consensus
        read = getTestRead(false);
        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = getTestRead(false);
        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = getTestRead(false);
        addConsensusReadAttribute(read, 2, 1, UmiReadType.SINGLE);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(3, combinedStats.readCounts().TotalReads);
        assertEquals(1, combinedStats.readCounts().Duplicates);
        assertEquals(3, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(20, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(10, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
    }

    @Test
    public void testDualStrandSupplementaries()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        BamReader bamReader = new BamReader(new ChrBaseRegion(
                CHR_1, 1, 1000), mConfig, null, null, combinedStats, null);

        // supplementary with a duplicate
        SAMRecord read = getTestRead(true);
        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = getTestRead(true);
        read.setDuplicateReadFlag(true);

        bamReader.processRead(read);

        read = getTestRead(true);
        addConsensusReadAttribute(read, 2, 1, UmiReadType.DUAL);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(2, combinedStats.readCounts().TotalReads);
        assertEquals(1, combinedStats.readCounts().Duplicates);
        assertEquals(1, combinedStats.readCounts().DualStrand);
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(2, combinedStats.flagStats().passCount(FlagStatType.SUPPLEMENTARY));
        assertEquals(10, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(10, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
    }

    @Test
    public void testDualStrandSecondaryConsensus()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        BamReader bamReader = new BamReader(new ChrBaseRegion(
                CHR_1, 1, 1000), mConfig, null, null, combinedStats, null);

        // secondary with a duplicate
        SAMRecord read = getTestRead(false);
        read.setDuplicateReadFlag(true);
        read.setSecondaryAlignment(true);

        bamReader.processRead(read);

        read = getTestRead(false);
        read.setDuplicateReadFlag(true);
        read.setSecondaryAlignment(true);

        bamReader.processRead(read);

        read = getTestRead(false);
        addConsensusReadAttribute(read, 2, 1, UmiReadType.DUAL);
        read.setSecondaryAlignment(true);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(2, combinedStats.readCounts().TotalReads);
        assertEquals(1, combinedStats.readCounts().Duplicates);
        assertEquals(1, combinedStats.readCounts().DualStrand);
        assertEquals(2, combinedStats.flagStats().passCount(FlagStatType.TOTAL));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.SUPPLEMENTARY));
        assertEquals(2, combinedStats.flagStats().passCount(FlagStatType.SECONDARY));
        assertEquals(0, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(0, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
    }

    @Test
    public void testSingleStrandLowMappingQualityConsensus()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        BamReader bamReader = new BamReader(new ChrBaseRegion(
                CHR_1, 1, 1000), mConfig, null, null, combinedStats, null);

        SAMRecord read = getTestRead(false);
        read.setDuplicateReadFlag(true);
        read.setMappingQuality(0);

        bamReader.processRead(read);

        read = getTestRead(false);
        read.setDuplicateReadFlag(true);
        read.setMappingQuality(0);

        bamReader.processRead(read);

        read = getTestRead(false);
        addConsensusReadAttribute(read, 2, 1, UmiReadType.SINGLE);
        read.setMappingQuality(0);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(2, combinedStats.readCounts().TotalReads);
        assertEquals(1, combinedStats.readCounts().Duplicates);
        assertEquals(0, combinedStats.readCounts().DualStrand);
        assertEquals(2, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.SUPPLEMENTARY));
        assertEquals(0, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(0, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
        assertEquals(20, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.LOW_MAP_QUAL.ordinal()]);
    }

    @Test
    public void testSingleStrandLowBaseQualityConsensus()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        BamReader bamReader = new BamReader(new ChrBaseRegion(
                CHR_1, 1, 1000), mConfig, null, null, combinedStats, null);

        byte[] lowBaseQualities = new byte[10];
        Arrays.fill(lowBaseQualities, (byte) 5);

        SAMRecord read = getTestRead(false);
        read.setDuplicateReadFlag(true);
        read.setBaseQualities(lowBaseQualities);

        bamReader.processRead(read);

        read = getTestRead(false);
        read.setDuplicateReadFlag(true);
        read.setBaseQualities(lowBaseQualities);

        bamReader.processRead(read);

        read = getTestRead(false);
        addConsensusReadAttribute(read, 2, 1, UmiReadType.SINGLE);
        read.setBaseQualities(lowBaseQualities);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(2, combinedStats.readCounts().TotalReads);
        assertEquals(1, combinedStats.readCounts().Duplicates);
        assertEquals(0, combinedStats.readCounts().DualStrand);
        assertEquals(2, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.SUPPLEMENTARY));
        assertEquals(0, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(10, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
        assertEquals(10, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.LOW_BASE_QUAL.ordinal()]);
    }

    @Test
    public void testPartialOverlappingReadWithoutPreviousRegion()
    {
        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);
        BamReader bamReader = new BamReader(new ChrBaseRegion(
                CHR_1, 25, 50), mConfig, null, null, combinedStats, null);

        SAMRecord read = getTestRead(false);

        bamReader.processRead(read);

        bamReader.postSliceProcess();
        assertEquals(1, combinedStats.readCounts().TotalReads);
        assertEquals(0, combinedStats.readCounts().Duplicates);
        assertEquals(0, combinedStats.readCounts().DualStrand);
        assertEquals(1, combinedStats.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(0, combinedStats.flagStats().passCount(FlagStatType.SUPPLEMENTARY));
        assertEquals(5, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(0, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
        assertEquals(0, combinedStats.coverageMetrics().FilterTypeCounts[FilterType.LOW_BASE_QUAL.ordinal()]);
    }

    @Test
    public void testPartialOverlappingReadWithPreviousRegion()
    {
        CombinedStats combinedStatsWithoutPreviousRegion = new CombinedStats(mConfig.MaxCoverage);
        BamReader bamReaderWithoutPreviousRegion = new BamReader(
                new ChrBaseRegion(CHR_1, 25, 50),
                mConfig,
                null,
                null,
                combinedStatsWithoutPreviousRegion,
                new ChrBaseRegion(CHR_1, 1, 24)
                );

        SAMRecord read = getTestRead(false);

        bamReaderWithoutPreviousRegion.processRead(read);

        bamReaderWithoutPreviousRegion.postSliceProcess();
        assertEquals(0, combinedStatsWithoutPreviousRegion.readCounts().TotalReads);
        assertEquals(0, combinedStatsWithoutPreviousRegion.readCounts().Duplicates);
        assertEquals(0, combinedStatsWithoutPreviousRegion.readCounts().DualStrand);
        assertEquals(0, combinedStatsWithoutPreviousRegion.flagStats().passCount(FlagStatType.PRIMARY));
        assertEquals(0, combinedStatsWithoutPreviousRegion.flagStats().passCount(FlagStatType.DUPLICATE));
        assertEquals(0, combinedStatsWithoutPreviousRegion.flagStats().passCount(FlagStatType.PRIMARY_DUPLICATE));
        assertEquals(0, combinedStatsWithoutPreviousRegion.flagStats().passCount(FlagStatType.SUPPLEMENTARY));
        assertEquals(0, combinedStatsWithoutPreviousRegion.coverageMetrics().FilterTypeCounts[FilterType.UNFILTERED.ordinal()]);
        assertEquals(0, combinedStatsWithoutPreviousRegion.coverageMetrics().FilterTypeCounts[FilterType.DUPLICATE.ordinal()]);
        assertEquals(0, combinedStatsWithoutPreviousRegion.coverageMetrics().FilterTypeCounts[FilterType.LOW_BASE_QUAL.ordinal()]);
    }

    @NotNull
    private SAMRecord getTestRead(boolean isSupplementary)
    {
        return SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 20, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100,
                false, isSupplementary, null);
    }
}
