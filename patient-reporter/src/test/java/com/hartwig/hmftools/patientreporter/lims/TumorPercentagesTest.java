package com.hartwig.hmftools.patientreporter.lims;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class TumorPercentagesTest {

    private static final String RESOURCE_DIR = Resources.getResource("csv").getPath();
    private static final String TUMOR_PERCENTAGE_CSV = RESOURCE_DIR + File.separator + "tumor_percentages.csv";

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canLoadAndRetrieveTumorPercentages() throws IOException, EmptyFileException {
        final TumorPercentages tumorPercentages = TumorPercentages.loadFromCsv(TUMOR_PERCENTAGE_CSV);
        assertEquals(0.3, tumorPercentages.findTumorPercentageForSample("CPCT02020202T"), EPSILON);
        assertEquals(1, tumorPercentages.findTumorPercentageForSample("CPCT02030303T"), EPSILON);
        assertEquals(Double.NaN, tumorPercentages.findTumorPercentageForSample("CPCT02040404T"), EPSILON);
    }
}