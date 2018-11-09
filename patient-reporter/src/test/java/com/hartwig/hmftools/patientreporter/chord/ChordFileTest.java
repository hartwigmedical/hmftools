package com.hartwig.hmftools.patientreporter.chord;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;

import org.junit.Test;

public class ChordFileTest {

    private static final String CHORD_FILE = Resources.getResource("test_run/chord/CPCT11111111T_chord_prediction.txt").getPath();

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canReadTestChordPrediction() throws IOException {
        ChordAnalysis chordValues = ChordFileReader.read(CHORD_FILE);

        assertEquals(0.023, chordValues.hrdValue(), EPSILON);
    }
}