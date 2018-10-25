package com.hartwig.hmftools.patientreporter.chordclassifier;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class ChordFileTest {
    private static final String CHORD_FILE = Resources.getResource("chordclassifier/chord_prediction.txt").getPath();
    private static final double EPSILON = 1.0e-10;

    @Test
    public void canReadTestChordPrediction() throws IOException {
        List<ChordAnalysis> chordValues = ChordFile.loadChordFile(CHORD_FILE);

        assertEquals(1, chordValues.size());
        assertEquals(0.023, chordValues.iterator().next().hrdValue(), EPSILON);
    }
}