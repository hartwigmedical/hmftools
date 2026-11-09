package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.junit.Test;

public class VirusReportFileTest
{
    @Test
    public void testWriteReadRoundtrip() throws IOException
    {
        List<VirusReport> reports = List.of(
                new VirusReport(
                        "Human gammaherpesvirus 8", "Human gammaherpesvirus 8", ReportQcStatus.AMBIGUOUS_VIRUS, VirusType.HHV8,
                        2, 91.614, 19.8877, 2.6517, ClonalityStatus.CLONAL, true),
                new VirusReport(
                        "JC polyomavirus", "JC polyomavirus", ReportQcStatus.NO_ABNORMALITIES, null,
                        0, 83.5478, 14.9107, null, ClonalityStatus.PURPLE_FAIL_NO_TUMOR, false));

        File file = File.createTempFile("virus_report", ".tsv");
        file.deleteOnExit();
        VirusReportFile.write(file.getAbsolutePath(), reports);

        assertEquals(reports, VirusReportFile.read(file.getAbsolutePath()));
    }
}
