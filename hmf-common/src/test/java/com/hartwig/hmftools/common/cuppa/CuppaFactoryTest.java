package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class CuppaFactoryTest {

    private static final String CUPPA_DATA_CSV = Resources.getResource("cuppa/sample.cup.data.csv").getPath();

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canReadFromTestFile() throws IOException {
        List<CuppaEntry> entries = CuppaDataFile.read(CUPPA_DATA_CSV);

        CuppaData cuppa = CuppaFactory.create(entries);

        assertEquals("Esophagus/Stomach", cuppa.predictedCancerType());
        assertEquals(0.000132, cuppa.bestPredictionLikelihood(), EPSILON);
        assertEquals(10, cuppa.simpleDups32To200B());
        assertEquals(8, cuppa.maxComplexSize());
        assertEquals(5, cuppa.LINECount());
        assertEquals(0, cuppa.telomericSGLs());
    }

    @Test
    public void doNotCrashOnMissingEntries() {
        assertNotNull(CuppaFactory.create(Lists.newArrayList()));
    }
}