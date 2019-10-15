package com.hartwig.hmftools.common.variant.structural.linx;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class LinxDriverFileTest {
    private static final String LINX_DRIVERS_TSV = Resources.getResource("variant/structural/linx/sample.linx.drivers.tsv").getPath();

    @Test
    public void canReadLinxDriversTsv() throws IOException {
        List<LinxDriver> linxDriverTsv = LinxDriverFile.read(LINX_DRIVERS_TSV);
        assertEquals(2, linxDriverTsv.size());

        LinxDriver homozygousDisruption1 = linxDriverTsv.get(0);
        assertEquals(-50, homozygousDisruption1.clusterId());
        assertEquals("BRCA1", homozygousDisruption1.gene());
        assertEquals("LOH_CHR", homozygousDisruption1.eventType());

        LinxDriver homozygousDisruption2 = linxDriverTsv.get(1);
        assertEquals(270, homozygousDisruption2.clusterId());
        assertEquals("ATM", homozygousDisruption2.gene());
        assertEquals("HOM_DEL_DISRUPTION", homozygousDisruption2.eventType());
    }

}