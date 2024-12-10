package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

import htsjdk.samtools.SAMRecord;

public class HlaTransformerTest extends RemapperTestBase
{
    HlaRecordAligner aligner;
    HlaTransformer transformer;

    @Before
    public void setup()
    {
        aligner = Mockito.mock(HlaRecordAligner.class);
        transformer = new HlaTransformer(aligner);
    }

    @Test
    public void isHlaAltReferenceTest()
    {
        Assert.assertTrue(HlaTransformer.isHlaAltReference("HLA-A*01:03"));
        Assert.assertTrue(HlaTransformer.isHlaAltReference("hla-whatever"));
        Assert.assertFalse(HlaTransformer.isHlaAltReference("hla whatever"));
        Assert.assertFalse(HlaTransformer.isHlaAltReference("chrUn_KI270590v1"));
    }

    @Test
    public void hasAltReference()
    {
        Assert.assertFalse(HlaTransformer.hasAltReference(records.get(0)));
        Assert.assertFalse(HlaTransformer.hasAltReference(records.get(1)));
        Assert.assertTrue(HlaTransformer.hasAltReference(records.get(7)));
        Assert.assertTrue(HlaTransformer.hasAltReference(records.get(8)));
        Assert.assertTrue(HlaTransformer.hasAltReference(records.get(11)));
    }

    @Test
    public void mateHasAltReference()
    {
        Assert.assertFalse(HlaTransformer.mateHasAltReference(records.get(0)));
        Assert.assertTrue(HlaTransformer.mateHasAltReference(records.get(4)));
        Assert.assertTrue(HlaTransformer.mateHasAltReference(records.get(5)));
        Assert.assertFalse(HlaTransformer.mateHasAltReference(records.get(7)));
        Assert.assertTrue(HlaTransformer.mateHasAltReference(records.get(8)));
        Assert.assertTrue(HlaTransformer.mateHasAltReference(records.get(10)));
    }

    @Test
    public void hasSomeHlaReference()
    {
        Assert.assertFalse(HlaTransformer.hasSomeHlaReference(records.get(0)));
        Assert.assertTrue(HlaTransformer.hasSomeHlaReference(records.get(4)));
        Assert.assertTrue(HlaTransformer.hasSomeHlaReference(records.get(7)));
        Assert.assertTrue(HlaTransformer.hasSomeHlaReference(records.get(11)));
    }

    @Test
    public void aNonHlaRecordShouldPassThrough() {
        List<SAMRecord> processed = transformer.process(records.get(0));
        Assert.assertEquals(1, processed.size());
        Assert.assertEquals(records.get(0), processed.get(0));
        checkThatAlignerHasNotBeenCalled();
        checkThatNoRecordsRemainUnprocessed();
    }

    @Test
    public void hlaSupplementaryRecordsShouldBeIgnored() {
        Assert.assertTrue(HlaTransformer.mateHasAltReference(records.get(5)));
        Assert.assertTrue(records.get(5).isSecondaryOrSupplementary()); // sanity
        List<SAMRecord> processed = transformer.process(records.get(5));
        Assert.assertEquals(0, processed.size());
        checkThatAlignerHasNotBeenCalled();
        checkThatNoRecordsRemainUnprocessed();
    }

    @Test
    public void hlaPrimaryRecordsShouldBeHeldUntilAMatchingRecordIsProcessed()
    {
        // Records 8 and 9 are a pair and are both HLA, as are records 10 and 11.
        List<SAMRecord> mockedAlignmentResult_8_9= List.of(records.get(8), records.get(9));
        Mockito.when(aligner.alignPair(pair(9,8))).thenReturn(mockedAlignmentResult_8_9);
        List<SAMRecord> mockedAlignmentResult_10_11= List.of(records.get(10), records.get(11));
        Mockito.when(aligner.alignPair(pair(11,10))).thenReturn(mockedAlignmentResult_10_11);

        List<SAMRecord> processed = transformer.process(records.get(8));
        Assert.assertEquals(0, processed.size());
        checkThatAlignerHasNotBeenCalled();

        processed = transformer.process(records.get(11));
        Assert.assertEquals(0, processed.size());
        checkThatAlignerHasNotBeenCalled();

        processed = transformer.process(records.get(9));
        Assert.assertEquals(mockedAlignmentResult_8_9, processed);
        Mockito.verify(aligner).alignPair(pair(9,8));

        processed = transformer.process(records.get(10));
        Assert.assertEquals(mockedAlignmentResult_10_11, processed);
        Mockito.verify(aligner).alignPair(pair(11,10));

//        Mockito.verify(aligner, Mockito.times(0)).alignRecord(Mockito.any());
    }

    @Test
    public void unmatchedHlaRecordsCanBeProcessed()
    {
        transformer.process(records.get(8));
        transformer.process(records.get(11));

//        Mockito.when(aligner.alignRecord(Mockito.any())).thenReturn(List.of(records.get(0)));
//        List<SAMRecord> returned = transformer.unmatchedRecords();
//        Assert.assertEquals(2, returned.size());
//         We don't know the order in which the unmatched records are returned.
//        Assert.assertTrue(returned.contains(records.get(8)));
//        Assert.assertTrue(returned.contains(records.get(11)));
    }

    private RecordPair pair(int leftIndex, int rightIndex)
    {
        return new RecordPair(records.get(leftIndex), records.get(rightIndex));
    }

    private void checkThatAlignerHasNotBeenCalled()
    {
        Mockito.verify(aligner, Mockito.times(0)).alignPair(Mockito.any());
//        Mockito.verify(aligner, Mockito.times(0)).alignRecord(Mockito.any());
    }

    private void checkThatNoRecordsRemainUnprocessed()
    {
        Assert.assertEquals(0, transformer.unmatchedRecords().size());
    }
}
