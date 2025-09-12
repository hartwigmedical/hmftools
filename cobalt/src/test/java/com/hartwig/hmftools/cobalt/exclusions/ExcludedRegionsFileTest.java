package com.hartwig.hmftools.cobalt.exclusions;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

public class ExcludedRegionsFileTest
{
    public static String line(String chromosome, int start, int end, String reason)
    {
        return chromosome + "\t" + start + "\t" + end + "\t" + reason;
    }

    public static void writeExcludedRegionsFile(File file, List<String> lines) throws IOException
    {
        List<String> allLines = new ArrayList<>();
        allLines.add("chromosome\tstart\tend\treason");
        allLines.addAll(lines);
        FileUtils.writeLines(file, allLines);
    }

    @Test
    public void load() throws IOException
    {
        File tempDir = FileUtils.getTempDirectory();
        File regionsFile = new File(tempDir, "regions.tsv");
        regionsFile.deleteOnExit();
        List<String> lines = new ArrayList<>();
        lines.add(line("15", 22482842, 22482962, "IG_PSEUDO_REGION"));
        lines.add(line("16", 31973576, 31973696, "IG_PSEUDO_REGION"));
        lines.add(line("16", 32063475, 32063595, "IG_PSEUDO_REGION"));
        writeExcludedRegionsFile(regionsFile, lines);

        BufferedReader reader = new BufferedReader(new FileReader(regionsFile));
        List<ChrBaseRegion> loaded = new ExcludedRegionsFile(reader).regions();
        reader.close();
        assertEquals(3, loaded.size());
        assertEquals("15", loaded.get(0).Chromosome);
        assertEquals(22482842, loaded.get(0).start());
        assertEquals(22482962, loaded.get(0).end());
        assertEquals("16", loaded.get(1).Chromosome);
        assertEquals(31973576, loaded.get(1).start());
        assertEquals(31973696, loaded.get(1).end());
        assertEquals("16", loaded.get(2).Chromosome);
        assertEquals(32063475, loaded.get(2).start());
        assertEquals(32063595, loaded.get(2).end());
    }
}
