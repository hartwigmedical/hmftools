package com.hartwig.hmftools.amber.blacklist;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class AmberBlacklistFileTest
{
    @Test
    public void readWriteTest() throws IOException
    {
        List<AmberBlacklistPoint> statistics = new ArrayList<>();
        statistics.add(new AmberBlacklistPoint("chr1", 100, 200, 0.11));
        statistics.add(new AmberBlacklistPoint("chr1", 200, 100, 0.12));
        statistics.add(new AmberBlacklistPoint("chr2", 100, 10, 0.15));
        List<AmberBlacklistPoint> roundTripStatistics = writeThenRead(statistics);
        Assert.assertEquals(statistics, roundTripStatistics);
    }

    @Test
    public void doubleFormattingTest() throws IOException
    {
        List<AmberBlacklistPoint> statistics = new ArrayList<>();
        statistics.add(new AmberBlacklistPoint("chr1", 100, 200, 0.111111111));
        statistics.add(new AmberBlacklistPoint("chr1", 200, 100, 0.1211111));
        List<AmberBlacklistPoint> roundTripStatistics = writeThenRead(statistics);
        Assert.assertEquals(0.11, roundTripStatistics.get(0).meanVaf(), 0.000000001);
    }

    private List<AmberBlacklistPoint> writeThenRead(List<AmberBlacklistPoint> statistics) throws IOException
    {
        File tempDir = Files.createTempDirectory("amber").toFile();
        File destination = new File(tempDir, "statistics.tsv");
        if(destination.exists())
        {
            Files.delete(destination.toPath());
        }
        AmberBlacklistFile.writeToFile(destination, statistics);
        return AmberBlacklistFile.readFromFile(destination);
    }
}
