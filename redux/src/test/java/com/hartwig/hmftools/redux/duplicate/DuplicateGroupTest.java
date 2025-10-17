package com.hartwig.hmftools.redux.duplicate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createFragmentCoords;
import static com.hartwig.hmftools.redux.TestUtils.createTestConfig;
import static com.hartwig.hmftools.redux.TestUtils.setBaseQualities;

import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.duplicate.DuplicateGroup;
import com.hartwig.hmftools.redux.duplicate.DuplicateGroupBuilder;
import com.hartwig.hmftools.redux.duplicate.FragmentCoords;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroupTest
{
    @Test
    public void testNonConsensusPrimaryRead()
    {
        ReduxConfig config = createTestConfig();

        DuplicateGroupBuilder duplicateGroupBuilder = new DuplicateGroupBuilder(config);

        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();

        int readPos = 100;
        int matePos = 200;

        SAMRecord read1 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, matePos, false, false, null, true, TEST_READ_CIGAR);

        setBaseQualities(read1, 25);

        FragmentCoords fragmentCoords = createFragmentCoords(read1);

        SAMRecord read2 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, matePos, false, false, null, true, TEST_READ_CIGAR);

        setBaseQualities(read2, 25);

        SAMRecord read3 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, matePos, false, false, null, true, TEST_READ_CIGAR);

        setBaseQualities(read3, 11);

        DuplicateGroup duplicateGroup = new DuplicateGroup(List.of(read1, read2, read3), fragmentCoords);
        duplicateGroups.add(duplicateGroup);

        duplicateGroupBuilder.processDuplicateGroups(duplicateGroups, Collections.emptyList(), true);

        assertTrue(duplicateGroup.isPrimaryRead(read1));
    }
}
