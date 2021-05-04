package com.hartwig.hmftools.common.flagstat;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class FlagstatFileTest {

    private static final String FLAGSTAT_FILE = Resources.getResource("flagstat/example.flagstat").getPath();
    private static final String MALFORMED_FILE = Resources.getResource("flagstat/malformed.flagstat").getPath();

    private static final double EPSILON = 1E-10;

    @Test
    public void canReadFlagstatFile() throws IOException {
        Flagstat flagstat = FlagstatFile.loadFromFile(FLAGSTAT_FILE);

        assertEquals(0.1, flagstat.duplicateProportion(), EPSILON);
        assertEquals(0.8, flagstat.mappedProportion(), EPSILON);
    }

    @Test (expected = IOException.class)
    public void crashesOnMalformed() throws IOException {
        FlagstatFile.loadFromFile(MALFORMED_FILE);
    }
}