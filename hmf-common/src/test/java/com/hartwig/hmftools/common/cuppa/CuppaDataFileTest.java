package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class CuppaDataFileTest {

    private static final String CUPPA_DATA_CSV = Resources.getResource("cuppa/sample.cup.data.csv").getPath();

    @Test
    public void canLoadTestFile() throws IOException {
        List<CuppaEntry> entries = CuppaDataFile.read(CUPPA_DATA_CSV);
        assertEquals(29, entries.size());
    }
}