package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.junit.jupiter.api.Test;

import htsjdk.samtools.SAMRecord;

public class ReadCachePairedReadJitterAdaptorTest
{
    @Test
    public void testHandlesReadCachePossibleHoldingBackOfReadsForPerformanceReasons()
    {
        ReadCachePairedReadJitterAdaptor readCache =  new ReadCachePairedReadJitterAdaptor(new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, true, ILLUMINA));
        String readBases = "A".repeat(143);
        String cigar = "143M";

        // read1 and read2 need to be popped off together
        SAMRecord read1 = createSamRecord("READ_001", CHR_1, 850, readBases, cigar, CHR_1, 1047, false, false, null, true, cigar);
        SAMRecord read2 = createSamRecord("READ_002", CHR_1, 851, readBases, cigar, CHR_1, 1047, false, false, null, true, cigar);
        SAMRecord read3 = createSamRecord("READ_003", CHR_1, 1000, readBases, cigar, CHR_1, 919, true, false, null, false, cigar);
        SAMRecord read4 = createSamRecord("READ_004", CHR_1, 1098, readBases, cigar, CHR_1, 914, true, false, null, false, cigar);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3, read4);

        final Function<FragmentCoordReads, Stream<SAMRecord>> fragmentCoordReadsToReads = (final FragmentCoordReads fragmentCoordReads) -> Stream.concat(
                fragmentCoordReads.DuplicateGroups.stream().flatMap(x -> x.reads().stream()),
                fragmentCoordReads.SingleReads.stream().map(ReadInfo::read)
        );

        Set<String> readNamesOfInterest = Sets.newHashSet("READ_001", "READ_002");
        int totalPopsOfInterest = 0;
        for(SAMRecord read : reads)
        {
            readCache.processRead(read);
            FragmentCoordReads fragmentCoordReads = readCache.popReads();
            if(fragmentCoordReads == null)
                continue;

            int popsOfInterest = (int) fragmentCoordReadsToReads.apply(fragmentCoordReads).filter(x -> readNamesOfInterest.contains(x.getReadName())).count();
            totalPopsOfInterest += popsOfInterest;

            assertTrue(popsOfInterest == 0 || popsOfInterest == 2);
        }

        FragmentCoordReads fragmentCoordReads = readCache.evictAll();
        if(fragmentCoordReads != null)
        {
            int popsOfInterest = (int) fragmentCoordReadsToReads.apply(fragmentCoordReads).filter(x -> readNamesOfInterest.contains(x.getReadName())).count();
            totalPopsOfInterest += popsOfInterest;

            assertTrue(popsOfInterest == 0 || popsOfInterest == 2);
        }

        assertEquals(2, totalPopsOfInterest);
    }
}
