package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.RegionPerfTracker.BIN_SIZE;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

public class RegionPerfTrackerTest
{
    private static final long ONE_SECOND = 1_000_000_000L;

    @Test
    public void testPositionsInSameBinAggregate() throws IOException
    {
        RegionPerfTracker tracker = new RegionPerfTracker();
        tracker.add("chr1", 10, ONE_SECOND, 4);
        tracker.add("chr1", BIN_SIZE - 1, 2 * ONE_SECOND, 6);

        assertEquals(1, tracker.regionCount());

        Map<String, String[]> rows = writeAndRead(tracker);
        String[] row = rows.get("chr1:0");

        assertEquals("0", row[1]);
        assertEquals(String.valueOf(BIN_SIZE - 1), row[2]);
        assertEquals("2", row[3]);
        assertEquals("10", row[4]);
        assertEquals("3.000", row[5]);
        assertEquals("2.000", row[6]);
    }

    @Test
    public void testPositionsInDifferentBinsSeparate() throws IOException
    {
        RegionPerfTracker tracker = new RegionPerfTracker();
        tracker.add("chr1", 10, ONE_SECOND, 2);
        tracker.add("chr1", BIN_SIZE, ONE_SECOND, 2);
        tracker.add("chr2", 10, ONE_SECOND, 2);

        assertEquals(3, tracker.regionCount());

        Map<String, String[]> rows = writeAndRead(tracker);
        assertEquals(String.valueOf(BIN_SIZE), rows.get("chr1:" + BIN_SIZE)[1]);
        assertEquals("0", rows.get("chr2:0")[1]);
    }

    @Test
    public void testUnmappedGroupsBucketSeparately() throws IOException
    {
        RegionPerfTracker tracker = new RegionPerfTracker();
        tracker.add(null, 0, ONE_SECOND, 2);
        tracker.add("chr1", 10, ONE_SECOND, 2);

        assertEquals(2, tracker.regionCount());
        assertEquals("2", writeAndRead(tracker).get("unmapped:0")[4]);
    }

    @Test
    public void testMergeCombinesTotalsAndTakesMax() throws IOException
    {
        RegionPerfTracker first = new RegionPerfTracker();
        first.add("chr1", 10, ONE_SECOND, 3);

        RegionPerfTracker second = new RegionPerfTracker();
        second.add("chr1", 20, 5 * ONE_SECOND, 7);
        second.add("chr3", 20, ONE_SECOND, 1);

        first.merge(second);

        assertEquals(2, first.regionCount());

        String[] row = writeAndRead(first).get("chr1:0");
        assertEquals("2", row[3]);
        assertEquals("10", row[4]);
        assertEquals("6.000", row[5]);
        assertEquals("5.000", row[6]);
    }

    @Test
    public void testMergeLeavesSourceUnchanged() throws IOException
    {
        RegionPerfTracker first = new RegionPerfTracker();
        first.add("chr1", 10, ONE_SECOND, 3);

        RegionPerfTracker second = new RegionPerfTracker();
        second.add("chr1", 10, ONE_SECOND, 3);

        first.merge(second);
        first.merge(second);

        assertEquals("3", writeAndRead(second).get("chr1:0")[4]);
        assertEquals("9", writeAndRead(first).get("chr1:0")[4]);
    }

    private static Map<String, String[]> writeAndRead(final RegionPerfTracker tracker) throws IOException
    {
        Path file = Files.createTempFile("region_perf", ".tsv");
        tracker.write(file.toString(), Double.MAX_VALUE);

        List<String> lines = Files.readAllLines(file);
        Files.deleteIfExists(file);

        Map<String, String[]> rows = new HashMap<>();
        for(int i = 1; i < lines.size(); ++i)
        {
            String[] values = lines.get(i).split("\t", -1);
            rows.put(values[0] + ":" + values[1], values);
        }

        return rows;
    }
}
