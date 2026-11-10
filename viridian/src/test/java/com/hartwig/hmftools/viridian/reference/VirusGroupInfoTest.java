package com.hartwig.hmftools.viridian.reference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Map;

import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Test;

public class VirusGroupInfoTest
{
    private static final String HEADER = "virus_group\treporting_type";

    @Test
    public void testLoadsValidResource() throws IOException
    {
        // The HHV-8 reporting type differs from its enum constant name.

        Map<VirusGroup, VirusGroupInfo> info = load(
                "Human gammaherpesvirus 4\tEBV",
                "Human gammaherpesvirus 8\tHHV-8",
                "BK polyomavirus\tnull");

        VirusGroup ebv = new VirusGroup("Human gammaherpesvirus 4");
        VirusGroup hhv8 = new VirusGroup("Human gammaherpesvirus 8");
        VirusGroup bk = new VirusGroup("BK polyomavirus");
        Map<VirusGroup, VirusGroupInfo> expected = Map.of(
                ebv, new VirusGroupInfo(ebv, VirusType.EBV),
                hhv8, new VirusGroupInfo(hhv8, VirusType.HHV8),
                bk, new VirusGroupInfo(bk, null));
        assertEquals(expected, info);
    }

    @Test
    public void testDuplicateGroupFails()
    {
        assertThrows(
                UserInputError.class, () -> load(
                        "BK polyomavirus\tnull",
                        "BK polyomavirus\tnull"));
    }

    @Test
    public void testUnknownReportingTypeFails()
    {
        assertThrows(RuntimeException.class, () -> load("BK polyomavirus\tWART"));
    }

    private static Map<VirusGroup, VirusGroupInfo> load(String... rows) throws IOException
    {
        return VirusGroupInfo.load(writeFile(HEADER, rows).toString());
    }

    private static Path writeFile(String header, String... rows) throws IOException
    {
        Path file = Files.createTempFile("virus_groups", ".tsv");
        file.toFile().deleteOnExit();
        StringBuilder content = new StringBuilder(header).append('\n');
        Arrays.stream(rows).forEach(row -> content.append(row).append('\n'));
        Files.writeString(file, content.toString());
        return file;
    }
}
