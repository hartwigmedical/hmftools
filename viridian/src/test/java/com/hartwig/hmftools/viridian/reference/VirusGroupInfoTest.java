package com.hartwig.hmftools.viridian.reference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Map;

import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class VirusGroupInfoTest
{
    private static final String HEADER = "virus_group\treporting_type";

    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

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

    private Map<VirusGroup, VirusGroupInfo> load(String... rows) throws IOException
    {
        File file = mTempDir.newFile("virus_groups.tsv");
        Files.writeString(file.toPath(), HEADER + "\n" + String.join("\n", rows) + "\n");
        return VirusGroupInfo.load(file.getPath());
    }
}
