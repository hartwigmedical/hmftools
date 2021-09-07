package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class CuppaFactoryTest {

    private static final String CUPPA_DATA_CSV = Resources.getResource("cuppa/sample.cup.data.csv").getPath();

    @Test
    public void canReadFromTestFile() throws IOException {
        List<CuppaEntry> entries = CuppaDataFile.read(CUPPA_DATA_CSV);

        CuppaData cuppa = CuppaFactory.build(Strings.EMPTY, entries);

        assertEquals(10, cuppa.simpleDups32To200B());
        assertEquals(8, cuppa.maxComplexSize());
        assertEquals(5, cuppa.LINECount());
        assertEquals(0, cuppa.telomericSGLs());
    }
}