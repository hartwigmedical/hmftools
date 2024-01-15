package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.read.ReadFilters.recordSoftClipsNearJunction;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.esvee.common.Junction;

import org.junit.Test;

public class ReadFiltersTest
{
    @Test
    public void testBasicFilters()
    {
        Read read = createSamRecord("READ_01", 20, REF_BASES.substring(15, 48), "5S20M1S");

        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 19, POS_STRAND)));
        assertTrue(recordSoftClipsNearJunction(read, new Junction(CHR_1, 18, NEG_STRAND)));
        assertTrue(recordSoftClipsNearJunction(read, new Junction(CHR_1, 22, NEG_STRAND)));
        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 17, NEG_STRAND)));

        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 39, NEG_STRAND)));
        assertTrue(recordSoftClipsNearJunction(read, new Junction(CHR_1, 39, POS_STRAND)));
        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 40, POS_STRAND))); // doesn't extend far enough
    }


    /* CHASHA FIXME
    @Test
    public void canFilterLowQualityRecords()
    {
        final Record record = createSAMRecord("AAAA");
        record.getBaseQuality()[0] = 8;
        record.getBaseQuality()[1] = 11;
        record.getBaseQuality()[2] = 14;
        record.getBaseQuality()[3] = 11;

        assertTrue(AlignmentFilters.isRecordAverageQualityAbove(record, 12)).isFalse();
        assertTrue(AlignmentFilters.isRecordAverageQualityAbove(record, 11)).isFalse();
        assertTrue(AlignmentFilters.isRecordAverageQualityAbove(record, 10)).isTrue();
    }

    @Test
    public void canFilterRecordsWithLowQualityInRegionOfInterestForwards()
    {
        final Record record = createSAMRecord("AAAAA");
        record.getBaseQuality()[0] = 1;
        record.getBaseQuality()[1] = 2;
        record.getBaseQuality()[2] = 3;
        record.getBaseQuality()[3] = 4;
        record.getBaseQuality()[4] = 5;

        final Junction junction = mock(Junction.class);
        when(junction.position()).thenReturn(3);
        when(junction.orientation()).thenReturn(Direction.FORWARDS);

        assertTrue(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 3)).isTrue();
        assertTrue(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 4)).isFalse();
        assertTrue(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 5)).isFalse();
    }

    @Test
    public void canFilterRecordsWithLowQualityInRegionOfInterestReverse()
    {
        final Record record = createSAMRecord("AAAAA");
        record.getBaseQuality()[0] = 1;
        record.getBaseQuality()[1] = 2;
        record.getBaseQuality()[2] = 3;
        record.getBaseQuality()[3] = 4;
        record.getBaseQuality()[4] = 5;

        final Junction junction = mock(Junction.class);
        when(junction.position()).thenReturn(3);
        when(junction.orientation()).thenReturn(Direction.REVERSE);

        assertTrue(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 1)).isTrue();
        assertTrue(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 2)).isFalse();
        assertTrue(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 3)).isFalse();
    }

    @Test
    public void canFilterRecordsThatMatchReferenceForwards()
    {
        // This test assumes READ_SOFT_CLIP_JUNCTION_TOLERANCE is 2.

        final Junction junction = mock(Junction.class);
        when(junction.position()).thenReturn(2);
        when(junction.orientation()).thenReturn(Direction.FORWARDS);

        final Record record = createSAMRecord("AAAAAA", 1);
        record.setCigar("6M");
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("3M3S");
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isTrue();

        record.setCigar("5M1S");
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("4S2M"); // We shouldn't see data that looks like this, but if we do this is the current behaviour
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();
    }

    @Test
    public void canFilterRecordsThatMatchReferenceReverse()
    {
        // This test assumes READ_SOFT_CLIP_JUNCTION_TOLERANCE is 2.

        final Junction junction = mock(Junction.class);
        when(junction.position()).thenReturn(5);
        when(junction.orientation()).thenReturn(Direction.REVERSE);

        final Record record = createSAMRecord("AAAAAAA", 1);
        record.setCigar("7M");
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("4S3M");
        record.setAlignmentStart(1 + 3);
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isTrue();

        record.setCigar("1S5M");
        record.setAlignmentStart(1 + 1);
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("2M5S"); // We shouldn't see data that looks like this, but if we do this is the current behaviour
        assertTrue(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();
    }
    */
}
