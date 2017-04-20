package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class LimsTest {

    private static final String LIMS_EXAMPLE_FILE =
            Resources.getResource("lims").getPath() + File.separator + "LimsExample.csv";
    private static final double EPSILON = 1.0e-10;

    @Test
    public void canLoadFromCsv() throws IOException, EmptyFileException {
        final LimsModel limsModel = Lims.buildModelFromCsv(LIMS_EXAMPLE_FILE);

        final Map<String, LimsBiopsyData> dataPerSample = limsModel.data();
        assertEquals(2, dataPerSample.size());

        final String sample2Name = "SAMPLE2";
        final LimsBiopsyData sample2 = dataPerSample.get(sample2Name);
        assertNotNull(sample2);
        assertNotNull(sample2.samplingDate());
        assertNotNull(sample2.arrivalDate());
        assertEquals(0.8, sample2.tumorPercentage(), EPSILON);
        assertEquals(sample2.tumorPercentage(), limsModel.findTumorPercentageForSample(sample2Name), EPSILON);

        final String sample3Name = "SAMPLE3";
        final LimsBiopsyData sample3 = dataPerSample.get(sample3Name);
        assertNotNull(sample3);
        assertNull(sample3.samplingDate());
        assertNotNull(sample3.arrivalDate());
        assertEquals(0.3, sample3.tumorPercentage(), EPSILON);
        assertEquals(sample3.tumorPercentage(), limsModel.findTumorPercentageForSample(sample3Name), EPSILON);
    }
}
