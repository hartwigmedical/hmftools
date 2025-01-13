package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.Arrays;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class HlaAlignmentTest extends RemapperTestBase
{
    @Test
    public void isUnmapped()
    {
        BwaMemAlignment bwa = bwa("181,-1,-1,-1,-1,-1,0,0,0,0,\"\",null,null,5,29943927,0");
        HlaAlignment alignment = new HlaAlignment(bwa);
        Assert.assertTrue(alignment.isUnmapped());

        BwaMemAlignment bwa2 = bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0");
        HlaAlignment alignment2 = new HlaAlignment(bwa2);
        Assert.assertFalse(alignment2.isUnmapped());
    }

    @Test
    public void createSAMRecord()
    {
        String readName = "A00624:8:HHKYHDSXX:1:1446:18213:29684";
        String nukes = "GGCCAGGGTCTCACACCCTCCAGAGCATGTACGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATAACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCCTGGACCGCCGCGGACACGGC";
        String qualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        RawFastaData data = new RawFastaData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719");
        HlaAlignment alignment = new HlaAlignment(bwa);

        BwaMemAlignment bwa2 = bwa("145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719");
        HlaAlignment mate = new HlaAlignment(bwa2);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        Assert.assertEquals(readName, sam.getReadName());
        Assert.assertEquals(719, sam.getInferredInsertSize());
        Assert.assertFalse(sam.getReadNegativeStrandFlag());
        Assert.assertArrayEquals(nukes.getBytes(), sam.getReadBases());
        Assert.assertArrayEquals(qualities.getBytes(), sam.getBaseQualities());
        Assert.assertEquals(mate.getRefId(), sam.getReferenceIndex());
        Assert.assertEquals(31355730, sam.getAlignmentStart());
        Assert.assertEquals(mate.Position_1Based, sam.getMateAlignmentStart());
        Assert.assertEquals(97, sam.getFlags());
        Assert.assertEquals(60, sam.getMappingQuality());
        Assert.assertEquals("151M", sam.getCigarString());
    }

    @Test
    public void createSAMRecordReverseStrand()
    {
        String readName = "A00624:8:HHKYHDSXX:1:1446:18213:29684";
        String nukes = "TTGTTCTCTGCCTCACACTCAGTGTGTTTGGGGCTCTGATTCCAGCACTTCTGAGTCACTTTACCTCCACTCAGATCAGGAGCAGAAGTCCCTGTTCCCCGCTCAGAGACTCGAACTTTCCAATGAATAGGAGATTATCCCAGGTGCCTGC";
        String qualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF";
        RawFastaData data = new RawFastaData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719");
        HlaAlignment alignment = new HlaAlignment(bwa);

        BwaMemAlignment bwa2 = bwa("97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719");
        HlaAlignment mate = new HlaAlignment(bwa2);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        Assert.assertEquals(readName, sam.getReadName());
        Assert.assertEquals(-719, sam.getInferredInsertSize());
        Assert.assertTrue(sam.getReadNegativeStrandFlag());
        Assert.assertArrayEquals(Nucleotides.reverseComplementBases(nukes.getBytes()), sam.getReadBases());
        Assert.assertArrayEquals(Arrays.reverseArray(qualities.getBytes()), sam.getBaseQualities());
        Assert.assertEquals(mate.getRefId(), sam.getReferenceIndex());
        Assert.assertEquals(31356298, sam.getAlignmentStart());
        Assert.assertEquals(mate.Position_1Based, sam.getMateAlignmentStart());
        Assert.assertEquals(145, sam.getFlags());
        Assert.assertEquals(60, sam.getMappingQuality());
        Assert.assertEquals("151M", sam.getCigarString());
    }

    @Test
    public void createSAMRecordUnmapped()
    {
        String readName = "A00624:8:HHKYHDSXX:4:2543:5737:19288";
        String nukes = "TTTTTTTTTTTACTTTTTTTGTTTATGTTTGGTGGGTGTAACCTATGTTTCATTTTTTGTTTATGTTTTTTTAAACCATGGTGTTTTCCTCGCCTTCTGGACCTGTGGTTGGTTTTTTTTTTTTTGTTGTTCCATTTTTTCCCCCCTCGGT";
        String qualities = ":FFF,FFFFF:,,FFF,F,,,,,F,F,F,::::FFF,F:F,:,F,F,F,F,F,F,FFF,FF,,F:FFF,F,,FF:,,,F,F,,,FF::,,,,F,,:,,FF,,,F,,:,,F,,:FF,F:FF,FF,,,F,F,,F,,,FF,,F,,FF,:,:,:F";
        RawFastaData data = new RawFastaData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("181,-1,-1,-1,-1,-1,0,0,0,0,,null,null,5,29943927,0");
        HlaAlignment alignment = new HlaAlignment(bwa);
        BwaMemAlignment bwa2 = bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0");
        HlaAlignment mate = new HlaAlignment(bwa2);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        Assert.assertEquals(readName, sam.getReadName());
        Assert.assertTrue(sam.getReadNegativeStrandFlag());
        Assert.assertArrayEquals(Nucleotides.reverseComplementBases(nukes.getBytes()), sam.getReadBases());
        Assert.assertArrayEquals(Arrays.reverseArray(qualities.getBytes()), sam.getBaseQualities());
        Assert.assertEquals(mate.getRefId(), sam.getReferenceIndex());
        Assert.assertEquals(mate.Position_1Based, sam.getAlignmentStart());
        Assert.assertEquals(mate.Position_1Based, sam.getMateAlignmentStart());
        Assert.assertEquals(181, sam.getFlags());
        Assert.assertEquals(0, sam.getMappingQuality());
        Assert.assertEquals("*", sam.getCigarString());
        Assert.assertEquals(0, sam.getInferredInsertSize());
    }

    @Test
    public void createSAMRecordMateUnmapped()
    {
        String readName = "A00624:8:HHKYHDSXX:4:2543:5737:19288";
        // A00624:8:HHKYHDSXX:4:2543:5737:19288	121	chr6	29943928	60	151M	=	29943928	FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:151	MQ:i:0	AS:i:151	XS:i:73
        String nukes = "GGACCAGAAGTCGCTGTTCCCTTCTCAGGGAATAGAAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCTCTTCCCCATCCCGGGTGTCCTGTCCATTCTCAAGATGGCCACATGCGTGCTGGTGGAGTGTCCCA";
        String qualities = "FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        RawFastaData data = new RawFastaData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0");
        HlaAlignment alignment = new HlaAlignment(bwa);

        BwaMemAlignment bwaMate = bwa("181,-1,-1,-1,-1,-1,0,0,0,0,,null,null,5,29943927,0");
        HlaAlignment mate = new HlaAlignment(bwaMate);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        Assert.assertEquals(readName, sam.getReadName());
        Assert.assertTrue(sam.getReadNegativeStrandFlag());
        Assert.assertArrayEquals(Nucleotides.reverseComplementBases(nukes.getBytes()), sam.getReadBases());
        Assert.assertArrayEquals(Arrays.reverseArray(qualities.getBytes()), sam.getBaseQualities());
        Assert.assertEquals(bwa.getRefId(), sam.getReferenceIndex().intValue());
        Assert.assertEquals(bwa.getRefStart() + 1, sam.getAlignmentStart()); // bwa record is 0-based
        Assert.assertEquals(bwa.getRefStart() + 1, sam.getMateAlignmentStart()); // NOT the value (0) from the raw bwa mate alignment
        Assert.assertEquals(121, sam.getFlags());
        Assert.assertEquals(60, sam.getMappingQuality());
        Assert.assertEquals("151M", sam.getCigarString());
        Assert.assertEquals(bwa.getRefId(), sam.getMateReferenceIndex().intValue());
        Assert.assertEquals("*", sam.getAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE));
        Assert.assertEquals(0, sam.getAttribute(SamRecordUtils.MATE_QUALITY_ATTRIBUTE));
        Assert.assertEquals(0, sam.getInferredInsertSize());
    }


}
