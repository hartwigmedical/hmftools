package com.hartwig.hmftools.cobalt.count;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class DepthReadingsFileTest
{
    @Test
    public void writeTest() throws Exception
    {
        List<DepthReading> readings = new ArrayList<>();
        readings.add(dr("1", 1, 88.34, 0.45));
        readings.add(dr("1", 1001, 80.04, 0.35));
        readings.add(dr("1", 2001, 85.43, 0.439));
        readings.add(dr("2", 1, 35.94, 0.45));
        readings.add(dr("2", 1001, 38.03, 0.50));
        readings.add(dr("3", 1, 182.82, 0.59));
        readings.add(dr("3", 1001, 188.10, 0.60));

        File destination = Files.createTempDirectory("readings").toFile();
        destination.deleteOnExit();
        File output = new File(destination, "depthReadings.tsv.gz");
        DepthReadingsFile.write(output.getAbsolutePath(), readings);

        List<DepthReading> recovered = new DepthReadingsFile(output.getAbsolutePath()).read();
        Assert.assertEquals(readings, recovered);
    }

    private DepthReading dr(String chr, int pos, double depth, double gc)
    {
        return new DepthReading(chr, pos, depth, gc);
    }
}
