package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_6;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.TEST_READ_ID;

public class RemapperTest {

    @Test
    public void hasAltReferenceTest() {
        SAMRecord noAltNoMate = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_6, 1234, "", TEST_READ_CIGAR, false, false, null);
        Assert.assertFalse(Remapper.hasAltReference(noAltNoMate));

        SAMRecord noAltMateNoAlt = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_6, 1234, "", TEST_READ_CIGAR, CHR_6, 100, false,
                false, null, false, TEST_READ_CIGAR);
        Assert.assertFalse(Remapper.hasAltReference(noAltMateNoAlt));

//        SAMRecord altNoMate = SamRecordTestUtils.createSamRecordUnpaired(
//                TEST_READ_ID, "HLA-A*01:11N", 1234, "", TEST_READ_CIGAR, false, false, null);
//        Assert.assertTrue(Remapper.hasAltReference(noAltNoMate));

    }

    private void readTestFile() {

    }
}
