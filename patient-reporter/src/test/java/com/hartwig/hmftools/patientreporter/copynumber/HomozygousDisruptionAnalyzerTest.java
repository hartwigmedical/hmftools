package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriverFile;

import org.junit.Test;

public class HomozygousDisruptionAnalyzerTest {
    private static final String LINX_DRIVERS_TSV = Resources.getResource("test_run/linx/sample.linx.drivers.tsv").getPath();

    @Test
    public void canAnnotateDelDisruptions() throws IOException {
        List<LinxDriver> homozygousDisruptionTsv = LinxDriverFile.read(LINX_DRIVERS_TSV);
        assertEquals(2, homozygousDisruptionTsv.size());

        List<LinxDriver> delDisruptions = HomozygousDisruptionAnalyzer.extractDelDisruptions(homozygousDisruptionTsv);
        assertEquals(1, delDisruptions.size());

        LinxDriver homozygousDisruption1 = homozygousDisruptionTsv.get(1);
        assertEquals(270, homozygousDisruption1.clusterId());
        assertEquals("ATM", homozygousDisruption1.gene());
        assertEquals("HOM_DEL_DISRUPTION", homozygousDisruption1.eventType());
    }

}