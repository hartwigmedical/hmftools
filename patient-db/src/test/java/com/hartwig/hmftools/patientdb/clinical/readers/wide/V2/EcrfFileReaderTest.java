package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.junit.Test;

public class EcrfFileReaderTest
{

    @Test
    public void testReadPatientDataDoesNotThrow() throws IOException
    {
        var result = EcrfFileReader.readPatientData("readers/wide/dummy_patient_data.csv");

        assertEquals(2, result.size());
    }

    @Test
    public void testReadBiopsyDataDoesNotThrow() throws IOException
    {
        var result = EcrfFileReader.readBiopsyData("readers/wide/dummy_biopsy_data.csv");

        assertEquals(2, result.size());
    }

    @Test
    public void testReadPrevTreatChemoDataDoesNotThrow() throws IOException
    {
        var result = EcrfFileReader.readPrevTreatChemoData("readers/wide/dummy_prev_treat_chemo_data.csv");

        assertEquals(1, result.size());
    }

    @Test
    public void testReadPrevTreatRadDoesNotThrow() throws IOException
    {
        var result = EcrfFileReader.readPrevTreatRadData("readers/wide/dummy_prev_treat_rad_data.csv");

        assertEquals(1, result.size());
    }

    @Test
    public void testReadTumorMeasureDoesNotThrow() throws IOException
    {
        var result = EcrfFileReader.readTumorMeasureData("readers/wide/dummy_tumor_measure_data.csv");

        assertEquals(8, result.size());
    }

    @Test
    public void testTreatChemoAvlDataDoesNotThrow() throws IOException
    {
        var result = EcrfFileReader.readTreatChemoAvlData("readers/wide/dummy_treat_chemo_avl_data.csv");

        assertEquals(1, result.size());
    }
}
