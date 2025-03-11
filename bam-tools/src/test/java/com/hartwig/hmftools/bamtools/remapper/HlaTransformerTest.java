package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.bamtools.remapper.HlaTransformer.isHlaAltReference;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

import htsjdk.samtools.SAMRecord;

public class HlaTransformerTest extends RemapperTestBase
{
    private HlaRecordAligner mAligner;
    private HlaTransformer mTransformer;

    @Before
    public void setup()
    {
        mAligner = Mockito.mock(HlaRecordAligner.class);
        mTransformer = new HlaTransformer(mAligner);
    }

    @Test
    public void isHlaAltReferenceTest()
    {
        assertTrue(isHlaAltReference("HLA-A*01:03"));
        assertTrue(isHlaAltReference("hla-whatever"));
        assertFalse(isHlaAltReference("hla whatever"));
        assertFalse(isHlaAltReference("chrUn_KI270590v1"));
    }

    @Test
    public void hasAltReference()
    {
        assertFalse(HlaTransformer.hasAltReference(records.get(0)));
        assertFalse(HlaTransformer.hasAltReference(records.get(1)));
        assertTrue(HlaTransformer.hasAltReference(records.get(7)));
        assertTrue(HlaTransformer.hasAltReference(records.get(8)));
        assertTrue(HlaTransformer.hasAltReference(records.get(11)));
    }

    @Test
    public void mateHasAltReference()
    {
        assertFalse(HlaTransformer.mateHasAltReference(records.get(0)));
        assertTrue(HlaTransformer.mateHasAltReference(records.get(4)));
        assertTrue(HlaTransformer.mateHasAltReference(records.get(5)));
        assertFalse(HlaTransformer.mateHasAltReference(records.get(7)));
        assertTrue(HlaTransformer.mateHasAltReference(records.get(8)));
        assertTrue(HlaTransformer.mateHasAltReference(records.get(10)));
    }

    @Test
    public void hasSomeHlaReference()
    {
        assertFalse(HlaTransformer.hasSomeHlaReference(records.get(0)));
        assertTrue(HlaTransformer.hasSomeHlaReference(records.get(4)));
        assertTrue(HlaTransformer.hasSomeHlaReference(records.get(7)));
        assertTrue(HlaTransformer.hasSomeHlaReference(records.get(11)));
    }

    @Test
    public void aNonHlaRecordShouldPassThrough()
    {
        List<SAMRecord> processed = mTransformer.process(records.get(0));
        Assert.assertEquals(1, processed.size());
        Assert.assertEquals(records.get(0), processed.get(0));
        checkThatAlignerHasNotBeenCalled();
        checkThatNoRecordsRemainUnprocessed();
    }

    @Test
    public void hlaSupplementaryRecordsShouldBeIgnored()
    {
        assertTrue(HlaTransformer.mateHasAltReference(records.get(5)));
        assertTrue(records.get(5).isSecondaryOrSupplementary()); // sanity
        List<SAMRecord> processed = mTransformer.process(records.get(5));
        Assert.assertEquals(0, processed.size());
        checkThatAlignerHasNotBeenCalled();
        checkThatNoRecordsRemainUnprocessed();
    }

    @Test
    public void hlaPrimaryRecordsShouldBeHeldUntilAMatchingRecordIsProcessed()
    {
        // Records 8 and 9 are a pair and are both HLA, as are records 10 and 11.
        List<SAMRecord> mockedAlignmentResult_8_9 = List.of(records.get(8), records.get(9));
        Mockito.when(mAligner.alignPair(pair(9, 8))).thenReturn(mockedAlignmentResult_8_9);
        List<SAMRecord> mockedAlignmentResult_10_11 = List.of(records.get(10), records.get(11));
        Mockito.when(mAligner.alignPair(pair(11, 10))).thenReturn(mockedAlignmentResult_10_11);

        List<SAMRecord> processed = mTransformer.process(records.get(8));
        Assert.assertEquals(0, processed.size());
        checkThatAlignerHasNotBeenCalled();

        processed = mTransformer.process(records.get(11));
        Assert.assertEquals(0, processed.size());
        checkThatAlignerHasNotBeenCalled();

        processed = mTransformer.process(records.get(9));
        Assert.assertEquals(mockedAlignmentResult_8_9, processed);
        Mockito.verify(mAligner).alignPair(pair(9, 8));

        processed = mTransformer.process(records.get(10));
        Assert.assertEquals(mockedAlignmentResult_10_11, processed);
        Mockito.verify(mAligner).alignPair(pair(11, 10));
    }

    @Test
    public void unmatchedHlaRecordsCanBeProcessed()
    {
        mTransformer.process(records.get(8));
        mTransformer.process(records.get(11));
    }

    private RecordPair pair(int leftIndex, int rightIndex)
    {
        return new RecordPair(records.get(leftIndex), records.get(rightIndex));
    }

    private void checkThatAlignerHasNotBeenCalled()
    {
        Mockito.verify(mAligner, Mockito.times(0)).alignPair(Mockito.any());
    }

    private void checkThatNoRecordsRemainUnprocessed()
    {
        Assert.assertEquals(0, mTransformer.unmatchedRecords().size());
    }
}
