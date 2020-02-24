package com.hartwig.hmftools.sage.read;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextDistanceTest {

    @Test
    public void testComplicated() {

        final String refBases = "ATTACAAT";
        final String readBases = "TTTTTAGTAGCCTTT";
        final String readCigar = "5S4M2I1M2D1M2S";
        final SAMRecord record = buildSamRecord(readCigar, readBases);

        ReadContextDistance victim = new ReadContextDistance(true, 1, 6, record, refBases.getBytes());
        assertDistance(1, "1M1X1M", victim);

        victim = new ReadContextDistance(true, 2, 6, record, refBases.getBytes());
        assertDistance(2, "1S1M1X2M", victim);

        victim = new ReadContextDistance(true, 3, 6, record, refBases.getBytes());
        assertDistance(3, "2S1M1X2M1I", victim);

        victim = new ReadContextDistance(true, 4, 6, record, refBases.getBytes());
        assertDistance(3, "3S1M1X2M2I", victim);

        victim = new ReadContextDistance(true, 5, 6, record, refBases.getBytes());
        assertDistance(3, "4S1M1X2M2I1M", victim);

        victim = new ReadContextDistance(true, 6, 6, record, refBases.getBytes());
        assertDistance(4, "5S1M1X2M2I1M2D1M", victim);
    }

    @Test
    public void testOnlyAltDifferent() {
        final String refBases  = "ATTACAAT";
        final String readBases = "ATTAGAAT";
        final String readCigar = "8M";

        final SAMRecord record = buildSamRecord(readCigar, readBases);
        final ReadContextDistance victim = new ReadContextDistance(true, 3, 4, record, refBases.getBytes());
        assertDistance(1, "3M1X3M", victim);
    }

    @Test
    public void testOnlyAltAndNeighbourDifferent() {
        final String refBases = "ATTACAAT";
        final String readBases = "ATTAGTAT";
        final String readCigar = "8M";

        final SAMRecord record = buildSamRecord(readCigar, readBases);
        final ReadContextDistance victim = new ReadContextDistance(true, 3, 4, record, refBases.getBytes());
        assertDistance(2, "3M2X2M", victim);
    }

    @Test
    public void testOverlapSoftClip() {
        final String refBases = "GATCA";
        final String readBases = "TTTTGGTCA";
        final String readCigar = "4S5M";

        final SAMRecord record = buildSamRecord(readCigar, readBases);
        final ReadContextDistance victim = new ReadContextDistance(true, 3, 5, record, refBases.getBytes());
        assertDistance(2, "2S1M1X3M", victim);
    }



    private void assertDistance(int expectedDistance, String expectedCigar, ReadContextDistance distance) {
        assertEquals(expectedCigar, distance.cigar());
        assertEquals(expectedDistance, distance.distance());
    }


    @NotNull
    static SAMRecord buildSamRecord(@NotNull final String cigar, @NotNull final String readString) {
        final StringBuilder qualityString = new StringBuilder();
        for (int i = 0; i < readString.length(); i++) {
            qualityString.append("A");
        }

        return buildSamRecord(cigar, readString, qualityString.toString());
    }

    @NotNull
     static SAMRecord buildSamRecord(@NotNull final String cigar, @NotNull final String readString, @NotNull final String qualities) {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(100);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        return record;
    }

}
