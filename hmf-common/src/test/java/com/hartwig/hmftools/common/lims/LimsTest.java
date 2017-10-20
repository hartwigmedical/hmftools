package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class LimsTest {

    private static final String LIMS_EXAMPLE_FILE = Resources.getResource("lims").getPath() + File.separator + "LimsExample.csv";
    private static final double EPSILON = 1.0e-10;
    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @Test
    public void canLoadFromCsv() throws IOException, EmptyFileException {
        final LimsModel limsModel = Lims.buildModelFromCsv(LIMS_EXAMPLE_FILE, dateFormatter);

        final Map<String, LimsData> dataPerSample = limsModel.data();
        assertEquals(4, dataPerSample.size());

        final String sample3RName = "SAMP01010003R";
        final LimsData sample3R = dataPerSample.get(sample3RName);
        assertNotNull(sample3R);
        assertNull(sample3R.samplingDate());
        assertNotNull(sample3R.arrivalDate());
        assertEquals("2016-01-01", sample3R.arrivalDate().toString());

        final String sample3TName = "SAMP01010003T";
        final LimsData sample3T = dataPerSample.get(sample3TName);
        final LocalDate sample3TSamplingDate = sample3T.samplingDate();
        assertNotNull(sample3T);
        assertNotNull(sample3TSamplingDate);
        assertNotNull(sample3T.arrivalDate());
        assertEquals("2016-01-02", sample3TSamplingDate.toString());
        assertEquals("2016-01-03", sample3T.arrivalDate().toString());
        assertTrue(sample3T.isTumor());
        final LimsTumorData sample3Tumor = (LimsTumorData) sample3T;
        assertNull(sample3Tumor.tumorPercentage());

        final String sample4Name = "SAMP01010004T";
        final LimsData sample4T = dataPerSample.get(sample4Name);
        final LocalDate sample4TSamplingDate = sample4T.samplingDate();
        assertNotNull(sample4T);
        assertNotNull(sample4TSamplingDate);
        assertNotNull(sample4T.arrivalDate());
        assertEquals("2016-01-04", sample4TSamplingDate.toString());
        assertEquals("2016-02-05", sample4T.arrivalDate().toString());
        assertTrue(sample4T.isTumor());
        final LimsTumorData sample4Tumor = (LimsTumorData) sample4T;
        final Double sample4TumorPercentage = sample4Tumor.tumorPercentage();
        assertNotNull(sample4TumorPercentage);
        assertEquals(0.8, sample4TumorPercentage, EPSILON);

        final String sample5Name = "SAMP01010005T";
        final LimsData sample5 = dataPerSample.get(sample5Name);
        assertNotNull(sample5);
        assertNull(sample5.samplingDate());
        assertNotNull(sample5.arrivalDate());
        assertEquals("2016-01-06", sample5.arrivalDate().toString());
        assertTrue(sample5.isTumor());
        final LimsTumorData sample5Tumor = (LimsTumorData) sample5;
        final Double sample5TumorPercentage = sample5Tumor.tumorPercentage();
        assertNotNull(sample5TumorPercentage);
        assertEquals(0.3, sample5TumorPercentage, EPSILON);
    }
}
