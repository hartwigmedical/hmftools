package com.hartwig.hmftools.patientreporter.chord;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ChordFileTest {

    private static final String CHORD_FILE = Resources.getResource("test_run/chord/CPCT11111111T_chord_prediction.txt").getPath();

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canReadTestChordPrediction() throws IOException {
        ChordAnalysis chordValues = ChordFile.fromFile(CHORD_FILE);

        assertEquals(0.023, chordValues.hrdValue(), EPSILON);
    }
}