package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Mockito;

import htsjdk.samtools.SAMRecord;

public class BwaHlaRecordAlignerTest extends RemapperTestBase
{

    private Aligner aligner;
    private BwaHlaRecordAligner bwaAligner;

    @Before
    public void before()
    {
        aligner = Mockito.mock(Aligner.class);
        bwaAligner = new BwaHlaRecordAligner(aligner);
    }

    @Test
    public void alignRecordTest()
    {
        List<BwaMemAlignment> alignments = List.of(
                bwa("16,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0"),
                bwa("2048,1,32916297,32916333,67,103,0,0,36,36,67S36M48S,36,null,-1,-1,0")
        );
        Mockito.when(aligner.alignSequence(Mockito.any())).thenReturn(alignments);
        SAMRecord record = records.get(13);
        List<SAMRecord> returned = bwaAligner.alignRecord(record);
        Assert.assertEquals(2, returned.size());
        Assert.assertEquals(HlaRecordAligner.createRemappedRecord(record, alignments.get(0)), returned.get(0));
        Assert.assertEquals(HlaRecordAligner.createRemappedRecord(record, alignments.get(1)), returned.get(1));

        ArgumentCaptor<byte[]> captor = ArgumentCaptor.forClass(byte[].class);
        Mockito.verify(aligner, Mockito.times(1)).alignSequence(captor.capture());
        Assert.assertArrayEquals(record.getReadBases(), captor.getValue());
    }

    @Test
    public void alignPairTest()
    {
        List<BwaMemAlignment> alignmentsForRecord12 = List.of(
                bwa("16,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,-1,-1,0")
        );
        List<BwaMemAlignment> alignmentsForRecord13 = List.of(
                bwa("16,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0"),
                bwa("2048,1,32916297,32916333,67,103,0,0,36,36,67S36M48S,36,null,-1,-1,0")
        );
        ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignments = ImmutablePair.of(alignmentsForRecord13, alignmentsForRecord12);
        Mockito.when(aligner.alignSequences(Mockito.any(), Mockito.any())).thenReturn(alignments);

        SAMRecord record12 = records.get(12);
        SAMRecord record13 = records.get(13);

        List<SAMRecord> returned = bwaAligner.alignPair(new RecordPair(record13, record12));
        Assert.assertEquals(3, returned.size());

        SAMRecord results0 = returned.get(0);
        // results0 comes from record13 realigned using the first element of alignmentsForRecord13
        // and with mate reference information set to match the single element of alignmentsForRecord12.
        final SAMRecord expected0 = HlaRecordAligner.createRemappedRecord(record13, alignmentsForRecord13.get(0));
        expected0.setMateReferenceIndex(5);
        expected0.setMateAlignmentStart(31354513);
        check(expected0, results0);

        SAMRecord results1 = returned.get(1);
        final SAMRecord expected1 = HlaRecordAligner.createRemappedRecord(record13, alignmentsForRecord13.get(1));
        expected1.setMateReferenceIndex(5);
        expected1.setMateAlignmentStart(31354513);
        check(expected1, results1);

        SAMRecord results2 = returned.get(2);
        // results2 comes from record12 realigned using the only element of alignmentsForRecord12
        // and with mate reference information set to match the first element of alignmentsForRecord13.
        final SAMRecord expected2 = HlaRecordAligner.createRemappedRecord(record12, alignmentsForRecord12.get(0));
        expected2.setMateReferenceIndex(5);
        expected2.setMateAlignmentStart(31354375);
        check(expected2, results2);

        ArgumentCaptor<byte[]> captor1 = ArgumentCaptor.forClass(byte[].class);
        ArgumentCaptor<byte[]> captor2 = ArgumentCaptor.forClass(byte[].class);
        Mockito.verify(aligner, Mockito.times(1)).alignSequences(captor1.capture(), captor2.capture());
        Assert.assertArrayEquals(record13.getReadBases(), captor1.getValue()); // Note that the pair is (record13, record12).
        Assert.assertArrayEquals(record12.getReadBases(), captor2.getValue());
    }
}
