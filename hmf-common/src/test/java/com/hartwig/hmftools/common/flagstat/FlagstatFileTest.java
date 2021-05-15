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
        Flagstat flagstat = FlagstatFile.read(FLAGSTAT_FILE);

        assertEquals(970, flagstat.uniqueReadCount());
        assertEquals(10, flagstat.secondaryCount());
        assertEquals(20, flagstat.supplementaryCount());

        assertEquals(0.1, flagstat.duplicateProportion(), EPSILON);
        assertEquals(0.8, flagstat.mappedProportion(), EPSILON);
        assertEquals(0.051546, flagstat.pairedInSequencingProportion(), EPSILON);
        assertEquals(0.927835, flagstat.properlyPairedProportion(), EPSILON);
        assertEquals(0.5, flagstat.withItselfAndMateMappedProportion(), EPSILON);
        assertEquals(0.04, flagstat.singletonProportion(), EPSILON);
    }

    @Test (expected = IOException.class)
    public void crashesOnMalformed() throws IOException {
        FlagstatFile.read(MALFORMED_FILE);
    }
}