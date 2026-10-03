package com.hartwig.hmftools.bamtools.remapper;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.utils.Arrays;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class HlaAlignmentTest extends RemapperTestBase
{
    @Test
    public void createSAMRecordFromTagData()
    {
        String readName = "A00624:8:HHKYHDSXX:1:1446:18213:29684";
        String nukes =
                "GGCCAGGGTCTCACACCCTCCAGAGCATGTACGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATAACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCCTGGACCGCCGCGGACACGGC";
        String qualities =
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        BamReadData data = new BamReadData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("97,5,31355729,31355880,0,151,59,0,151,19,151M,151,null,5,31356297,719");
        AltAlignment altAlignment = new AltAlignment("chrX", 12345, Orientation.FORWARD, "151M", 120);
        HlaAlignment alignment = new HlaAlignment(bwa, altAlignment);

        BwaMemAlignment bwa2 = bwa("145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719");
        HlaAlignment mate = new HlaAlignment(bwa2);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        assertEquals(readName, sam.getReadName());
        assertEquals(0, sam.getInferredInsertSize());
        Assert.assertFalse(sam.getReadNegativeStrandFlag());
        assertArrayEquals(nukes.getBytes(), sam.getReadBases());
        assertArrayEquals(qualities.getBytes(), sam.getBaseQualities());
        assertEquals(22, (int) sam.getReferenceIndex());
        assertEquals(12345, sam.getAlignmentStart());
        assertEquals(mate.Position, sam.getMateAlignmentStart());
        assertEquals(97, sam.getFlags());
        assertEquals(59, sam.getMappingQuality());
        assertEquals("151M", sam.getCigarString());
    }

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
        String nukes =
                "GGCCAGGGTCTCACACCCTCCAGAGCATGTACGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATAACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCCTGGACCGCCGCGGACACGGC";
        String qualities =
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        BamReadData data = new BamReadData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719");
        HlaAlignment alignment = new HlaAlignment(bwa);

        BwaMemAlignment bwa2 = bwa("145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719");
        HlaAlignment mate = new HlaAlignment(bwa2);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        assertEquals(readName, sam.getReadName());
        assertEquals(719, sam.getInferredInsertSize());
        Assert.assertFalse(sam.getReadNegativeStrandFlag());
        assertArrayEquals(nukes.getBytes(), sam.getReadBases());
        assertArrayEquals(qualities.getBytes(), sam.getBaseQualities());
        assertEquals(5, (int) sam.getReferenceIndex());
        assertEquals(31355730, sam.getAlignmentStart());
        assertEquals(mate.Position, sam.getMateAlignmentStart());
        assertEquals(97, sam.getFlags());
        assertEquals(60, sam.getMappingQuality());
        assertEquals("151M", sam.getCigarString());
    }

    @Test
    public void createSAMRecordReverseStrand()
    {
        String readName = "A00624:8:HHKYHDSXX:1:1446:18213:29684";
        String nukes =
                "TTGTTCTCTGCCTCACACTCAGTGTGTTTGGGGCTCTGATTCCAGCACTTCTGAGTCACTTTACCTCCACTCAGATCAGGAGCAGAAGTCCCTGTTCCCCGCTCAGAGACTCGAACTTTCCAATGAATAGGAGATTATCCCAGGTGCCTGC";
        String qualities =
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF";
        BamReadData data = new BamReadData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719");
        HlaAlignment alignment = new HlaAlignment(bwa);

        BwaMemAlignment bwa2 = bwa("97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719");
        HlaAlignment mate = new HlaAlignment(bwa2);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        assertEquals(readName, sam.getReadName());
        assertEquals(-719, sam.getInferredInsertSize());
        Assert.assertTrue(sam.getReadNegativeStrandFlag());
        assertArrayEquals(Nucleotides.reverseComplementBases(nukes.getBytes()), sam.getReadBases());
        assertArrayEquals(Arrays.reverseArray(qualities.getBytes()), sam.getBaseQualities());
        assertEquals(mate.getRefId(), sam.getReferenceIndex());
        assertEquals(31356298, sam.getAlignmentStart());
        assertEquals(mate.Position, sam.getMateAlignmentStart());
        assertEquals(145, sam.getFlags());
        assertEquals(60, sam.getMappingQuality());
        assertEquals("151M", sam.getCigarString());
    }

    @Test
    public void createSAMRecordUnmapped()
    {
        String readName = "A00624:8:HHKYHDSXX:4:2543:5737:19288";
        String nukes =
                "TTTTTTTTTTTACTTTTTTTGTTTATGTTTGGTGGGTGTAACCTATGTTTCATTTTTTGTTTATGTTTTTTTAAACCATGGTGTTTTCCTCGCCTTCTGGACCTGTGGTTGGTTTTTTTTTTTTTGTTGTTCCATTTTTTCCCCCCTCGGT";
        String qualities =
                ":FFF,FFFFF:,,FFF,F,,,,,F,F,F,::::FFF,F:F,:,F,F,F,F,F,F,FFF,FF,,F:FFF,F,,FF:,,,F,F,,,FF::,,,,F,,:,,FF,,,F,,:,,F,,:FF,F:FF,FF,,,F,F,,F,,,FF,,F,,FF,:,:,:F";
        BamReadData data = new BamReadData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("181,-1,-1,-1,-1,-1,0,0,0,0,,null,null,5,29943927,0");
        HlaAlignment alignment = new HlaAlignment(bwa);
        BwaMemAlignment bwa2 = bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0");
        HlaAlignment mate = new HlaAlignment(bwa2);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        assertEquals(readName, sam.getReadName());
        Assert.assertTrue(sam.getReadNegativeStrandFlag());
        assertArrayEquals(Nucleotides.reverseComplementBases(nukes.getBytes()), sam.getReadBases());
        assertArrayEquals(Arrays.reverseArray(qualities.getBytes()), sam.getBaseQualities());
        assertEquals(mate.getRefId(), sam.getReferenceIndex());
        assertEquals(mate.Position, sam.getAlignmentStart());
        assertEquals(mate.Position, sam.getMateAlignmentStart());
        assertEquals(181, sam.getFlags());
        assertEquals(0, sam.getMappingQuality());
        assertEquals("*", sam.getCigarString());
        assertEquals(0, sam.getInferredInsertSize());
    }

    @Test
    public void createSAMRecordMateUnmapped()
    {
        String readName = "A00624:8:HHKYHDSXX:4:2543:5737:19288";
        // A00624:8:HHKYHDSXX:4:2543:5737:19288	121	chr6	29943928	60	151M	=	29943928	FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:151	MQ:i:0	AS:i:151	XS:i:73
        String nukes =
                "GGACCAGAAGTCGCTGTTCCCTTCTCAGGGAATAGAAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCTCTTCCCCATCCCGGGTGTCCTGTCCATTCTCAAGATGGCCACATGCGTGCTGGTGGAGTGTCCCA";
        String qualities =
                "FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        BamReadData data = new BamReadData(readName, nukes.getBytes(), qualities.getBytes());
        BwaMemAlignment bwa = bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0");
        HlaAlignment alignment = new HlaAlignment(bwa);

        BwaMemAlignment bwaMate = bwa("181,-1,-1,-1,-1,-1,0,0,0,0,,null,null,5,29943927,0");
        HlaAlignment mate = new HlaAlignment(bwaMate);

        SAMRecord sam = alignment.createSamRecord(samFileHeader(), data, mate);
        assertEquals(readName, sam.getReadName());
        Assert.assertTrue(sam.getReadNegativeStrandFlag());
        assertArrayEquals(Nucleotides.reverseComplementBases(nukes.getBytes()), sam.getReadBases());
        assertArrayEquals(Arrays.reverseArray(qualities.getBytes()), sam.getBaseQualities());
        assertEquals(bwa.getRefId(), sam.getReferenceIndex().intValue());
        assertEquals(bwa.getRefStart() + 1, sam.getAlignmentStart()); // bwa record is 0-based
        assertEquals(bwa.getRefStart() + 1, sam.getMateAlignmentStart()); // NOT the value (0) from the raw bwa mate alignment
        assertEquals(121, sam.getFlags());
        assertEquals(60, sam.getMappingQuality());
        assertEquals("151M", sam.getCigarString());
        assertEquals(bwa.getRefId(), sam.getMateReferenceIndex().intValue());
        Assert.assertNull(sam.getAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE));
        assertEquals(0, sam.getAttribute(SamRecordUtils.MATE_QUALITY_ATTRIBUTE));
        assertEquals(0, sam.getInferredInsertSize());
    }

    @Test
    public void isUnmappedTest()
    {
        HlaAlignment unmapped = new HlaAlignment(bwa("181,-1,-1,-1,-1,-1,0,0,0,0,\"\",null,null,5,29943927,0"));
        Assert.assertTrue(unmapped.isUnmapped());

        HlaAlignment mapped = new HlaAlignment(bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0"));
        Assert.assertFalse(mapped.isUnmapped());
    }
}
