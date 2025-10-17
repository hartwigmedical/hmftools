package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.redux.ReduxConstants.UNMAP_MIN_HIGH_DEPTH;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.mappability.UnmappingRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.duplicate.FragmentCoords;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.unmap.ReadUnmapper;
import com.hartwig.hmftools.redux.unmap.UnmapRegionState;
import com.hartwig.hmftools.redux.write.BamWriter;

import htsjdk.samtools.SAMRecord;

public final class TestUtils
{
    public static final String TEST_READ_BASES = MockRefGenome.generateRandomBases(100);
    public static final String TEST_READ_ID = "READ_01";
    public static final String TEST_READ_CIGAR = "100M";

    public static final String REF_BASES_A = "AAAAAAAAAA";
    public static final String REF_BASES_C = "CCCCCCCCCC";
    public static final String REF_BASES_G = "GGGGGGGGGG";
    public static final String REF_BASES_T = "TTTTTTTTTT";

    public static final String CHR_4 = "4";
    public static final String CHR_5 = "5";

    public static final String REF_BASES_RANDOM = generateRandomBases(10);

    public static final String REF_BASES = "X" + REF_BASES_RANDOM + REF_BASES_A + REF_BASES_C + REF_BASES_G + REF_BASES_T
            + REF_BASES_A + REF_BASES_C + REF_BASES_G + REF_BASES_T + REF_BASES_RANDOM;

    public static final String REF_BASES_REPEAT_40 = REF_BASES.repeat(40);

    public static final int DEFAULT_QUAL = SamRecordTestUtils.DEFAULT_BASE_QUAL;

    public static final ReadIdGenerator READ_ID_GEN = new ReadIdGenerator();

    public static final ReadUnmapper READ_UNMAPPER_DISABLED = new ReadUnmapper(Collections.emptyMap());

    public static ReduxConfig createTestConfig()
    {
        return new ReduxConfig(
                new MockRefGenome(), false, false, false, READ_UNMAPPER_DISABLED);
    }

    public static PartitionReader createPartitionRead(final ReduxConfig config, final BamWriter writer)
    {
        PartitionReader partitionReader = new PartitionReader(config, null);
        partitionReader.setBamWriter(writer);
        return partitionReader;
    }

    public static void setBaseQualities(final SAMRecord read, int value)
    {
        for(int i = 0; i < read.getBaseQualities().length; ++i)
            read.getBaseQualities()[i] = (byte)value;
    }

    public static void setSecondInPair(final SAMRecord read)
    {
        read.setFirstOfPairFlag(false);
        read.setSecondOfPairFlag(true);
    }

    public static FragmentCoords createFragmentCoords(final SAMRecord read)
    {
        return FragmentCoords.fromRead(read, false);
    }

    public static ConsensusReadInfo createConsensusRead(
            final ConsensusReads consensusReads, final List<SAMRecord> reads, final String umiId)
    {
        FragmentCoords fragmentCoords = FragmentCoords.fromRead(reads.get(0), false);
        return consensusReads.createConsensusRead(reads, fragmentCoords, umiId);
    }

    public static SAMRecord createUnpairedRecord(final String readName, final String chromosome, final int readStart, int readEnd,
            boolean isReversed)
    {
        int readLength = readEnd - readStart + 1;
        String readBases = "A".repeat(readLength);
        String cigar = readLength + "M";
        return createSamRecordUnpaired(readName, chromosome, readStart, readBases, cigar, isReversed, false, null);
    }

    // unmapping test state
    public static final Map<String,List<UnmappingRegion>> CHR_LOCATION_MAP;
    public static final ReadUnmapper READ_UNMAPPER;

    static
    {
        CHR_LOCATION_MAP = Maps.newHashMap();
        CHR_LOCATION_MAP.put(CHR_1, Lists.newArrayList(new UnmappingRegion(500, 700, 0)));
        CHR_LOCATION_MAP.put(CHR_2, Lists.newArrayList());
        CHR_LOCATION_MAP.put(CHR_3, Lists.newArrayList(new UnmappingRegion(500, 700, UNMAP_MIN_HIGH_DEPTH)));

        CHR_LOCATION_MAP.put(CHR_4, Lists.newArrayList(
                new UnmappingRegion(1000, 2000, UNMAP_MIN_HIGH_DEPTH),
                new UnmappingRegion(3000, 4000, UNMAP_MIN_HIGH_DEPTH),
                new UnmappingRegion(5000, 6000, UNMAP_MIN_HIGH_DEPTH),
                new UnmappingRegion(7000, 8000, UNMAP_MIN_HIGH_DEPTH)));
        READ_UNMAPPER = new ReadUnmapper(CHR_LOCATION_MAP);
    }

    public static boolean checkTransformRead(final SAMRecord read, final String chromosome)
    {
        UnmapRegionState regionState = new UnmapRegionState(new ChrBaseRegion(
                chromosome, 1, 1000000), CHR_LOCATION_MAP.get(chromosome));

        return READ_UNMAPPER.checkTransformRead(read, regionState);
    }
}
