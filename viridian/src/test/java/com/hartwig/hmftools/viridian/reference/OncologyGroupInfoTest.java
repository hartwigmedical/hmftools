package com.hartwig.hmftools.viridian.reference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Map;

import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Test;

public class OncologyGroupInfoTest
{
    private static final String HEADER = "oncology_group\treporting_type\tdriver_likelihood";

    @Test
    public void testLoadsValidResource() throws IOException
    {
        Map<OncologyGroup, OncologyGroupInfo> info = load(
                "Human gammaherpesvirus 4\tEBV\tHIGH",
                "Human papillomavirus type 6\tHPV\tLOW",
                "BK polyomavirus\tnull\tnull");

        assertEquals(3, info.size());

        OncologyGroupInfo ebv = info.get(new OncologyGroup("Human gammaherpesvirus 4"));
        assertTrue(ebv.isReportable());
        assertEquals(VirusType.EBV, ebv.reportingType());
        assertEquals(VirusLikelihoodType.HIGH, ebv.driverLikelihood());

        OncologyGroupInfo bk = info.get(new OncologyGroup("BK polyomavirus"));
        assertFalse(bk.isReportable());
        assertNull(bk.reportingType());
        assertNull(bk.driverLikelihood());
    }

    @Test
    public void testHhv8ReportingTypeParses() throws IOException
    {
        Map<OncologyGroup, OncologyGroupInfo> info = load("Human gammaherpesvirus 8\tHHV-8\tHIGH");
        assertEquals(VirusType.HHV8, info.get(new OncologyGroup("Human gammaherpesvirus 8")).reportingType());
    }

    @Test
    public void testMissingColumnFails() throws IOException
    {
        Path file = writeFile("oncology_group\treporting_type", "BK polyomavirus\tnull");
        assertThrows(RuntimeException.class, () -> OncologyGroupInfo.load(file.toString()));
    }

    @Test
    public void testDuplicateGroupFails()
    {
        assertThrows(UserInputError.class, () -> load(
                "BK polyomavirus\tnull\tnull",
                "BK polyomavirus\tnull\tnull"));
    }

    @Test
    public void testReportingTypeWithoutDriverLikelihoodFails()
    {
        assertThrows(RuntimeException.class, () -> load("Human gammaherpesvirus 4\tEBV\tnull"));
    }

    @Test
    public void testDriverLikelihoodWithoutReportingTypeFails()
    {
        assertThrows(RuntimeException.class, () -> load("BK polyomavirus\tnull\tHIGH"));
    }

    @Test
    public void testUnknownReportingTypeFails()
    {
        assertThrows(RuntimeException.class, () -> load("BK polyomavirus\tWART\tHIGH"));
    }

    @Test
    public void testUnknownDriverLikelihoodFails()
    {
        assertThrows(RuntimeException.class, () -> load("Human gammaherpesvirus 4\tEBV\tMAYBE"));
    }

    @Test
    public void testUnknownDriverLikelihoodRejectsVirusLikelihoodUnknown()
    {
        // VirusLikelihoodType has UNKNOWN, but it is not a valid policy value.
        assertThrows(UserInputError.class, () -> load("Human gammaherpesvirus 4\tEBV\tUNKNOWN"));
    }

    private static Map<OncologyGroup, OncologyGroupInfo> load(String... rows) throws IOException
    {
        return OncologyGroupInfo.load(writeFile(HEADER, rows).toString());
    }

    private static Path writeFile(String header, String... rows) throws IOException
    {
        Path file = Files.createTempFile("oncology_group_info", ".tsv");
        file.toFile().deleteOnExit();
        StringBuilder content = new StringBuilder(header).append('\n');
        Arrays.stream(rows).forEach(row -> content.append(row).append('\n'));
        Files.writeString(file, content.toString());
        return file;
    }
}
