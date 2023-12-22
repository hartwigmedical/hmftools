package com.hartwig.hmftools.svassembly;

import static com.hartwig.hmftools.svassembly.TestUtils.createSAMRecord;

import static org.assertj.core.api.Assertions.assertThat;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import com.hartwig.hmftools.common.sv.Direction;
import com.hartwig.hmftools.svassembly.assembly.AlignmentFilters;
import com.hartwig.hmftools.svassembly.models.Record;

import org.junit.Test;

public class AlignmentFiltersTest
{
    @Test
    public void canFilterLowQualityRecords()
    {
        final Record record = createSAMRecord("AAAA");
        record.getBaseQuality()[0] = 8;
        record.getBaseQuality()[1] = 11;
        record.getBaseQuality()[2] = 14;
        record.getBaseQuality()[3] = 11;

        assertThat(AlignmentFilters.isRecordAverageQualityAbove(record, 12)).isFalse();
        assertThat(AlignmentFilters.isRecordAverageQualityAbove(record, 11)).isFalse();
        assertThat(AlignmentFilters.isRecordAverageQualityAbove(record, 10)).isTrue();
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

        assertThat(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 3)).isTrue();
        assertThat(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 4)).isFalse();
        assertThat(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 5)).isFalse();
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

        assertThat(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 1)).isTrue();
        assertThat(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 2)).isFalse();
        assertThat(AlignmentFilters.isRecordAverageQualityPastJunctionAbove(record, junction, 3)).isFalse();
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
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("3M3S");
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isTrue();

        record.setCigar("5M1S");
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("4S2M"); // We shouldn't see data that looks like this, but if we do this is the current behaviour
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();
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
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("4S3M");
        record.setAlignmentStart(1 + 3);
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isTrue();

        record.setCigar("1S5M");
        record.setAlignmentStart(1 + 1);
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();

        record.setCigar("2M5S"); // We shouldn't see data that looks like this, but if we do this is the current behaviour
        assertThat(AlignmentFilters.recordSoftClipsNearJunction(record, junction)).isFalse();
    }
}
