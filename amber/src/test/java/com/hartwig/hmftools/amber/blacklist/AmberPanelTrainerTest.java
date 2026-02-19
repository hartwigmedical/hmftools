package com.hartwig.hmftools.amber.blacklist;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class AmberPanelTrainerTest
{
    private File OutputFile;

    @Before
    public void setup() throws IOException
    {
        File TempDir = Files.createTempDirectory("amber").toFile();
        OutputFile = new File(TempDir, "output.tsv");
        if(OutputFile.exists())
        {
            Files.delete(OutputFile.toPath());
        }
    }

    @Test
    public void runTest() throws IOException
    {
        /*
        1	18099283	excluded - only in Sample3
        1	20026354	excluded - vaf ok
        2	208000822	average 0.32
        2	209099479	excluded - vaf ok
        3	128205519	average 0.7262
        3	128209482	average 0.7208
        4	87706617	excluded - vaf ok
        4	87706656	excluded - vaf ok
        5	131354359	excluded - vaf ok
        5	131354952	average 0.71
         */
        List<String> samples = List.of("Sample1", "Sample2", "Sample3", "Sample4");
        AmberPanelTrainer trainer = new AmberPanelTrainer(samples, testResourcesDir(), OutputFile);
        trainer.run();
        List<AmberBlacklistPoint> stats = AmberBlacklistFile.readFromFile(OutputFile);
        assertEquals(4, stats.size());
        assertEquals(new AmberBlacklistPoint("2", 208000822, 3, 0.32), stats.get(0));
        assertEquals(new AmberBlacklistPoint("3", 128205519, 3, 0.73), stats.get(1));
        assertEquals(new AmberBlacklistPoint("3", 128209482, 2, 0.72), stats.get(2));
        assertEquals(new AmberBlacklistPoint("5", 131354952, 3, 0.71), stats.get(3));
    }

    private File testResourcesDir()
    {
        return Path.of("src", "test", "resources", "blacklist").toFile();
    }
}
