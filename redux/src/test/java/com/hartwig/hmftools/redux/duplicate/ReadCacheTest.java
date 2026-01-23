package com.hartwig.hmftools.redux.duplicate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.ReduxConstants.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadCacheTest
{
    private static final String TEST_READ_BASES = "A".repeat(10);
    private static final String TEST_CIGAR = "10M";

    @Test
    public void testReadCacheAboveDynamicThreshold()
    {
        ReadCache readCache = new ReadCache(
                100, 50, false, new DuplicatesConfig(0),
                100, 100, 100);

        // create reads above the dynamic threshold but without standard purging occurring
        int readPos1 = 100;
        int readPos2 = 110;
        int readPos3 = 120;
        for(int i = 0; i < 40; ++i)
        {
            SAMRecord read = createSamRecord(
                    READ_ID_GEN.nextId(), CHR_1, readPos1, TEST_READ_BASES, TEST_CIGAR, CHR_1, 200, false, false,
                    null, true, TEST_CIGAR);

            readCache.processRead(read);

            read = createSamRecord(
                    READ_ID_GEN.nextId(), CHR_1, readPos2, TEST_READ_BASES, TEST_CIGAR, CHR_1, 200, false, false,
                    null, true, TEST_CIGAR);

            readCache.processRead(read);

            read = createSamRecord(
                    READ_ID_GEN.nextId(), CHR_1, readPos3, TEST_READ_BASES, TEST_CIGAR, CHR_1, 200, false, false,
                    null, true, TEST_CIGAR);

            readCache.processRead(read);
        }

        assertEquals(120, readCache.cachedReadCount());
        assertEquals(3, readCache.cachedFragCoordGroups());

        FragmentCoordReads fragmentCoordReads = readCache.popReads();
        assertNotNull(fragmentCoordReads);

        assertEquals(80, fragmentCoordReads.totalReadCount());
    }
}
