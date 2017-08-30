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
        final Double limsModelTumorPercentage4 = limsModel.findTumorPercentageForSample(sample4Name);
        assertNotNull(sample4TumorPercentage);
        assertNotNull(limsModelTumorPercentage4);
        assertEquals(0.8, sample4TumorPercentage, EPSILON);
        assertEquals(sample4TumorPercentage, limsModelTumorPercentage4, EPSILON);

        final String sample5Name = "SAMP01010005T";
        final LimsData sample5 = dataPerSample.get(sample5Name);
        assertNotNull(sample5);
        assertNull(sample5.samplingDate());
        assertNotNull(sample5.arrivalDate());
        assertEquals("2016-01-06", sample5.arrivalDate().toString());
        assertTrue(sample5.isTumor());
        final LimsTumorData sample5Tumor = (LimsTumorData) sample5;
        final Double sample5TumorPercentage = sample5Tumor.tumorPercentage();
        final Double limsModelTumorPercentage5 = limsModel.findTumorPercentageForSample(sample5Name);
        assertNotNull(sample5TumorPercentage);
        assertNotNull(limsModelTumorPercentage5);
        assertEquals(0.3, sample5TumorPercentage, EPSILON);
        assertEquals(sample5TumorPercentage, limsModelTumorPercentage5, EPSILON);
    }

    @Test
    public void canLoadJsonModel() throws IOException, EmptyFileException {
        final String limsJson = Resources.getResource("lims").getPath() + File.separator + "lims.json";
        final LimsJsonModel jsonLims = LimsJsonModel.readModelFromFile(limsJson);
        assertEquals(LocalDate.parse("2016-01-02"), jsonLims.samplingDateForSample("SAMP01010003R"));
        assertEquals(LocalDate.parse("2016-01-03"), jsonLims.arrivalDateForSample("SAMP01010003R"));
        assertNull(jsonLims.tumorPercentageForSample("SAMP01010003R"));
        assertEquals("CSB000000", jsonLims.bloodBarcodeForSample("SAMP01010003R"));
        assertEquals("CSB000000", jsonLims.barcodeForSample("SAMP01010003R"));

        assertEquals("CSB000000", jsonLims.bloodBarcodeForSample("SAMP01010003T"));
        assertEquals("FC0000001", jsonLims.barcodeForSample("SAMP01010003T"));
        final Double tumorPercentage = jsonLims.tumorPercentageForSample("SAMP01010003T");
        assertNotNull(tumorPercentage);
        assertEquals(30.0, tumorPercentage, 0.001);
        assertEquals(LocalDate.parse("2016-02-05"), jsonLims.arrivalDateForSample("SAMP01010003T"));
        assertEquals(LocalDate.parse("2016-01-04"), jsonLims.samplingDateForSample("SAMP01010003T"));
    }
}
