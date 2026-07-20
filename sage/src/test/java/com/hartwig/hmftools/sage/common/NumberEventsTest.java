package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class NumberEventsTest
{
    @Test
    public void testMissingNmForSplicedRead()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, "A".repeat(2_000));

        SAMRecord record = createSamRecordUnpaired(
                "read", CHR_1, 10, "AACAAAAC", "4M1000N4M", false, false, null);

        assertEquals(2, NumberEvents.calcAdjustedNumMutations(record, refGenome));

        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, 7);
        assertEquals(7, NumberEvents.calcAdjustedNumMutations(record, refGenome));
    }
}
