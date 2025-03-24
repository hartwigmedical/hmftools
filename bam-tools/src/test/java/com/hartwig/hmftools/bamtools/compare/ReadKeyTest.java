package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

public class ReadKeyTest
{
    @Test
    public void testCalcSupplementaryIndex()
    {
        // test that we can create keys that differentiate supplementary reads
        String readId = "read1";
        String readBases = "ACGT";
        String readCigar = "4M";

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 2000, readBases, readCigar, CHR_1, 2200,
                false, false, null);
        read.setMappingQuality(20);


        // this is not a supplementary, should have 0 as supplementary index
        ReadKey readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(0, readKey.SupplementaryIndex);

        // add supplementary attribute
        SupplementaryReadData suppData = new SupplementaryReadData(CHR_1, 2000, SUPP_POS_STRAND, readCigar, 0, 1);
        read.setAttribute(SAMTag.SA.name(), suppData.asSamTag());
        // read.setAttribute(SAMTag.SA.name(), "chr1,2000,+,151M,0,1;");

        // still since it is not supplementary, should have 0 as supplementary index
        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(0, readKey.SupplementaryIndex);

        // make it a supplementary
        read.setSupplementaryAlignmentFlag(true);

        // The SA attribute shows only the non supplementary read plus this read, so the supplementary index of this
        // read should be 1
        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(1, readKey.SupplementaryIndex);

        // set SA to show 2 supplementary reads, this time the other supplementary read is mapped to chr1:3000
        // which is > this read's genomic location, as a result the supplementary index is still 1
        suppData = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, readCigar, 0, 1);
        SupplementaryReadData suppData2 = new SupplementaryReadData(CHR_1, 3000, SUPP_POS_STRAND, readCigar, 0, 1);
        read.setAttribute(SAMTag.SA.name(), format("%s;%s", suppData.asSamTag(), suppData2.asSamTag()));

        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(1, readKey.SupplementaryIndex);

        // set SA to show 2 supplementary reads, this time the other supplementary read is mapped to chr1:1000
        // which is < this read's genomic location, as a result the supplementary index will be 2
        suppData2 = new SupplementaryReadData(CHR_1, 1500, SUPP_POS_STRAND, readCigar, 0, 1);
        read.setAttribute(SAMTag.SA.name(), format("%s;%s", suppData.asSamTag(), suppData2.asSamTag()));

        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(2, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, both are located genomically lower than this read
        // As a result the supplementary index will be 3
        SupplementaryReadData suppData3 = new SupplementaryReadData(CHR_1, 1800, SUPP_POS_STRAND, readCigar, 0, 1);

        read.setAttribute(SAMTag.SA.name(), format("%s;%s;%s", suppData.asSamTag(), suppData2.asSamTag(), suppData3.asSamTag()));

        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(3, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, all are located higher than this read
        suppData = new SupplementaryReadData(CHR_1, 3000, SUPP_POS_STRAND, readCigar, 0, 1);
        suppData2 = new SupplementaryReadData(CHR_1, 4000, SUPP_POS_STRAND, readCigar, 0, 1);
        suppData3 = new SupplementaryReadData(CHR_1, 5000, SUPP_POS_STRAND, readCigar, 0, 1);
        read.setAttribute(SAMTag.SA.name(), format("%s;%s;%s", suppData.asSamTag(), suppData2.asSamTag(), suppData3.asSamTag()));

        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(1, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, first one is located higher than this read, 2nd one is located lower
        suppData = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, readCigar, 0, 1);
        suppData2 = new SupplementaryReadData(CHR_1, 3000, SUPP_POS_STRAND, readCigar, 0, 1);
        suppData3 = new SupplementaryReadData(CHR_1, 1500, SUPP_POS_STRAND, readCigar, 0, 1);
        read.setAttribute(SAMTag.SA.name(), format("%s;%s;%s", suppData.asSamTag(), suppData2.asSamTag(), suppData3.asSamTag()));

        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(2, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, first one is located lower than this read, 2nd one is higher
        suppData = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, readCigar, 0, 1);
        suppData2 = new SupplementaryReadData(CHR_1, 1500, SUPP_POS_STRAND, readCigar, 0, 1);
        suppData3 = new SupplementaryReadData(CHR_1, 3000, SUPP_POS_STRAND, readCigar, 0, 1);
        read.setAttribute(SAMTag.SA.name(), format("%s;%s;%s", suppData.asSamTag(), suppData2.asSamTag(), suppData3.asSamTag()));

        readKey = ReadKey.from(read);
        Assert.assertEquals(readId, readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(2, readKey.SupplementaryIndex);
    }
}