package com.hartwig.hmftools.bamtools.remapper;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class HlaRecordAlignerTest extends RemapperTestBase
{
    @Test
    public void createRemappedRecordTest()
    {
        SAMRecord record = records.get(12);
        BwaMemAlignment alignment =
                new BwaMemAlignment(16, 5, 31354375, 31354419, 0, 44, 60, 1, 39, 20, "44M107S", "22C1", null, -1, -1, 0);
        SAMRecord remappedRecord = HlaRecordAligner.createRemappedRecord(record, alignment);

        Assert.assertEquals(record.getReadName(), remappedRecord.getReadName());
        Assert.assertArrayEquals(record.getReadBases(), remappedRecord.getReadBases());
        Assert.assertArrayEquals(record.getBaseQualities(), remappedRecord.getBaseQualities());
        Assert.assertEquals("chr6", remappedRecord.getReferenceName());
        Assert.assertEquals(alignment.getRefStart(), remappedRecord.getAlignmentStart());
        Assert.assertEquals(alignment.getMapQual(), remappedRecord.getMappingQuality());
    }
}
