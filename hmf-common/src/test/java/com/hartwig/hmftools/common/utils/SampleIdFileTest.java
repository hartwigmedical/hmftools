package com.hartwig.hmftools.common.utils;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

public class SampleIdFileTest
{
    @Test
    public void testSampleIdFileLoading()
    {
        List<String> sampleIds = loadDelimitedIdFile(loadSampleFile("sample_id_01.csv"), FLD_SAMPLE_ID, CSV_DELIM);
        assertEquals(3, sampleIds.size());
        assertEquals("SAMPLE_ID_01", sampleIds.get(0));
        assertEquals("SAMPLE_ID_03", sampleIds.get(2));

        sampleIds = loadDelimitedIdFile(loadSampleFile("sample_id_02.csv"), FLD_SAMPLE_ID, CSV_DELIM);
        assertEquals(3, sampleIds.size());
        assertEquals("SAMPLE_ID_01", sampleIds.get(0));
        assertEquals("SAMPLE_ID_03", sampleIds.get(2));

        sampleIds = loadDelimitedIdFile(loadSampleFile("sample_id_03.csv"), FLD_SAMPLE_ID, CSV_DELIM);
        assertEquals(3, sampleIds.size());
        assertEquals("SAMPLE_ID_01", sampleIds.get(0));
        assertEquals("SAMPLE_ID_03", sampleIds.get(2));

        // fails since headers are missing or incorrect
        sampleIds = loadDelimitedIdFile(loadSampleFile("sample_id_04.csv"), FLD_SAMPLE_ID, CSV_DELIM);
        assertTrue(sampleIds.isEmpty());

        sampleIds = loadDelimitedIdFile(loadSampleFile("sample_id_05.csv"), FLD_SAMPLE_ID, CSV_DELIM);
        assertTrue(sampleIds.isEmpty());
    }

    private static List<String> loadSampleFile(String sampleIdFile)
    {
        final InputStream inputStream = SampleIdFileTest.class.getResourceAsStream("/utils/" + sampleIdFile);
        return new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList());
    }
}
