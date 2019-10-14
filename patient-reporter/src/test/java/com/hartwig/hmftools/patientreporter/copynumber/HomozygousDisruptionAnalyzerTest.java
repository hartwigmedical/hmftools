package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.HomozygousDisruptionFile;

import org.junit.Test;

public class HomozygousDisruptionAnalyzerTest {
    private static final String LINX_DRIVERS_TSV = Resources.getResource("test_run/linx/sample.linx.drivers.tsv").getPath();

    @Test
    public void canAnnotateDelDisruptions() throws IOException {
        List<HomozygousDisruption> homozygousDisruptionTsv = HomozygousDisruptionFile.read(LINX_DRIVERS_TSV);
        assertEquals(2, homozygousDisruptionTsv.size());

        List<HomozygousDisruption> delDisruptions = HomozygousDisruptionAnalyzer.extractDelDisruptions(homozygousDisruptionTsv);
        assertEquals(1, delDisruptions.size());

        HomozygousDisruption homozygousDisruption1 = homozygousDisruptionTsv.get(1);
        assertEquals(270, homozygousDisruption1.clusterId());
        assertEquals("ATM", homozygousDisruption1.gene());
        assertEquals("HOM_DEL_DISRUPTION", homozygousDisruption1.eventType());
    }

}