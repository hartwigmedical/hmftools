package com.hartwig.hmftools.bamtools.remapper;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.remapper.BwaHlaRecordPairAligner.createSATag;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.test.ReadIdGenerator;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Mockito;

import htsjdk.samtools.SAMRecord;

public class BwaHlaRecordPairAlignerTest extends RemapperTestBase
{
    private final ReadIdGenerator mReadIdGen = new ReadIdGenerator();

    private PairAligner mAligner;
    private BwaHlaRecordPairAligner mBwaAligner;

    @Before
    public void before()
    {
        mAligner = Mockito.mock(PairAligner.class);
        mBwaAligner = new BwaHlaRecordPairAligner(mAligner, samFileHeader(), RefGenomeVersion.V38);
    }

    @Test
    public void alignPairTest()
    {
        List<BwaMemAlignment> alignmentsForRecord12 = List.of(
                bwa("145,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,5,31354375,-289")
        );
        List<BwaMemAlignment> alignmentsForRecord13 = List.of(
                bwa("97,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,5,31354513,289"),
                bwa("2161,1,32916486,32916522,67,103,0,0,36,36,67S36M48S,36,null,5,31354513,0"));
        ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignments =
                ImmutablePair.of(alignmentsForRecord13, alignmentsForRecord12);
        Mockito.when(mAligner.alignSequences(Mockito.any(), Mockito.any())).thenReturn(alignments);

        SAMRecord record12 = records.get(12);
        SAMRecord record13 = records.get(13);

        List<SAMRecord> returned = mBwaAligner.alignPair(new RecordPair(record13, record12));
        Assert.assertEquals(3, returned.size());

        ArgumentCaptor<byte[]> captor1 = ArgumentCaptor.forClass(byte[].class);
        ArgumentCaptor<byte[]> captor2 = ArgumentCaptor.forClass(byte[].class);
        Mockito.verify(mAligner, Mockito.times(1)).alignSequences(captor1.capture(), captor2.capture());
        Assert.assertArrayEquals(Nucleotides.reverseComplementBases(record13.getReadBases()), captor1.getValue()); // Note that the pair is (record13, record12).
        Assert.assertArrayEquals(record12.getReadBases(), captor2.getValue());
    }

    @Test
    public void addSupplementaryReadTagsTest()
    {
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int nmRead1 = 1;
        int supplementary1Pos = 1000;
        int nmSupplementary1 = 2;
        int supplementary2Pos = 2000;
        int nmSupplementary2 = 3;

        SAMRecord read1 = createSamRecord(nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                false,null, true, TEST_READ_CIGAR);
        read1.setAttribute(NUM_MUTATONS_ATTRIBUTE, nmRead1);

        SAMRecord supplementary1 = createSamRecord(read1.getReadName(), CHR_1, supplementary1Pos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, null, true, TEST_READ_CIGAR);
        supplementary1.setAttribute(NUM_MUTATONS_ATTRIBUTE, nmSupplementary1);

        SAMRecord supplementary2 = createSamRecord(
                read1.getReadName(), CHR_1, supplementary2Pos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, null, true, TEST_READ_CIGAR);
        supplementary2.setAttribute(NUM_MUTATONS_ATTRIBUTE, nmSupplementary2);

        List<SAMRecord> input = List.of(read1, supplementary1, supplementary2);
        BwaHlaRecordPairAligner.addSupplementaryReadTags(input);

        String saRead1 = read1.getAttribute(SUPPLEMENTARY_ATTRIBUTE).toString();
        String expectedSa1 = String.format("1,%d,+,100M,255,%d;1,%d,+,100M,255,%d;", supplementary1Pos, nmSupplementary1, supplementary2Pos, nmSupplementary2);
        Assert.assertEquals(expectedSa1, saRead1);

        String saSupplementary1 = supplementary1.getAttribute(SUPPLEMENTARY_ATTRIBUTE).toString();
        String expectedSaSupp1 = String.format("1,%d,+,100M,255,%d;1,%d,+,100M,255,%d;", readPos, nmRead1, supplementary2Pos, nmSupplementary2);
        Assert.assertEquals(expectedSaSupp1, saSupplementary1);

        String saSupplementary2 = supplementary2.getAttribute(SUPPLEMENTARY_ATTRIBUTE).toString();
        String expectedSaSupp2 = String.format("1,%d,+,100M,255,%d;1,%d,+,100M,255,%d;", readPos, nmRead1, supplementary1Pos, nmSupplementary1);
        Assert.assertEquals(expectedSaSupp2, saSupplementary2);
    }

    @Test
    public void createSATagTest()
    {
        // Records 5 and 8 are supplementary and primary alignments respectively for the same read.
        // Each contains an SA tag describing the other.
        // These can be used as expected values for the function that creates a tag from a record.
        // According to the specification, the SA tag of each should contain the cigar string of the other.
        // In practice, this is not always the case. In this example, the cigar in the SA tag of record8
        // is 42S32M77S but the cigar of record 5 is 42H32M77H.
        SAMRecord record5 = records.get(5);
        SAMRecord record8 = records.get(8);
        String saTag8 = record8.getAttribute(SUPPLEMENTARY_ATTRIBUTE).toString();
        String expected5 = saTag8.replace("42S32M77S", "42H32M77H");
        Assert.assertEquals(expected5, createSATag(record5));
        Assert.assertEquals(record5.getAttribute(SUPPLEMENTARY_ATTRIBUTE), createSATag(record8));
    }

    private String nextReadId(final String umiId)
    {
        return format("%s:%s", mReadIdGen.nextId(), umiId);
    }

}
