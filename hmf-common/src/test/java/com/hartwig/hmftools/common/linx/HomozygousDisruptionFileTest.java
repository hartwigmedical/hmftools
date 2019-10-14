package com.hartwig.hmftools.common.linx;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;

import org.junit.Test;

public class HomozygousDisruptionFileTest {
    private static final String LINX_DRIVERS_TSV = Resources.getResource("linx/sample.linx.drivers.tsv").getPath();

    @Test
    public void canReadLinxDriversTsv() throws IOException {
        List<HomozygousDisruption> homozygousDisruptionTsv = HomozygousDisruptionFile.read(LINX_DRIVERS_TSV);
        assertEquals(2, homozygousDisruptionTsv.size());

        HomozygousDisruption homozygousDisruption1 = homozygousDisruptionTsv.get(0);
        assertEquals(-50, homozygousDisruption1.clusterId());
        assertEquals("BRCA1", homozygousDisruption1.gene());
        assertEquals("LOH_CHR", homozygousDisruption1.eventType());

        HomozygousDisruption homozygousDisruption2 = homozygousDisruptionTsv.get(1);
        assertEquals(270, homozygousDisruption2.clusterId());
        assertEquals("ATM", homozygousDisruption2.gene());
        assertEquals("HOM_DEL_DISRUPTION", homozygousDisruption2.eventType());
    }
}