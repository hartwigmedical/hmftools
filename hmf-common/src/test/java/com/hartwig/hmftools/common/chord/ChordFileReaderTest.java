package com.hartwig.hmftools.common.chord;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ChordFileReaderTest {
    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final double EPSILON = 1.0E-10;

    private static final String SAMPLE = "CPCT11111111T";
    private static final double BRCA1 = 0.014;
    private static final double NONE = 0.894;
    private static final double BRCA2 = 0.050;
    private static final double HRD = 0.023;
    private static final double PREDICTED_RESPONSE = 0;

    @Test
    public void loadChordFile() throws IOException {
        String file = ChordFileReader.generateFilename(BASE_DIRECTORY, SAMPLE);

        ChordAnalysis chordAnalysis = ChordFileReader.read(file);

        assertEquals(BRCA1, chordAnalysis.BRCA1Value(), EPSILON);
        assertEquals(NONE, chordAnalysis.noneValue(), EPSILON);
        assertEquals(BRCA2, chordAnalysis.BRCA2Value(), EPSILON);
        assertEquals(HRD, chordAnalysis.hrdValue(), EPSILON);
        assertEquals(PREDICTED_RESPONSE, chordAnalysis.predictedResponseValue(), EPSILON);

        double brca1 = chordAnalysis.BRCA1Value();
        assertEquals(BRCA1, brca1, EPSILON);

        double none = chordAnalysis.noneValue();
        assertEquals(NONE, none, EPSILON);

        double brca2 = chordAnalysis.BRCA2Value();
        assertEquals(BRCA2, brca2, EPSILON);

        double hrd = chordAnalysis.hrdValue();
        assertEquals(HRD, hrd, EPSILON);

        double predicted_response = chordAnalysis.predictedResponseValue();
        assertEquals(PREDICTED_RESPONSE, predicted_response, EPSILON);

    }
}