package com.hartwig.hmftools.common.chord;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ChordFileReaderTest {
    private static final String FILE = Resources.getResource("chord/sample_chord_prediction.txt").getPath();
    private static final double EPSILON = 1.0E-10;

    private static final double NONE = 0.250;
    private static final double BRCA1 = 0.230;
    private static final double BRCA2 = 0.400;
    private static final double HRD = 0.630;
    private static final boolean PREDICTED_RESPONSE = true;

    @Test
    public void loadChordFile() throws IOException {
        ChordAnalysis chordAnalysis = ChordFileReader.read(FILE);

        assertEquals(NONE, chordAnalysis.noneValue(), EPSILON);
        assertEquals(BRCA1, chordAnalysis.BRCA1Value(), EPSILON);
        assertEquals(BRCA2, chordAnalysis.BRCA2Value(), EPSILON);
        assertEquals(HRD, chordAnalysis.hrdValue(), EPSILON);
        assertEquals(PREDICTED_RESPONSE, chordAnalysis.predictedResponseValue());
    }
}