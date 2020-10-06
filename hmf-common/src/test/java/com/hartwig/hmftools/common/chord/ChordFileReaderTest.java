package com.hartwig.hmftools.common.chord;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ChordFileReaderTest {

    private static final String V1_CHORD_FILE = Resources.getResource("chord/v1/sample_chord_prediction.txt").getPath();
    private static final String V2_CHORD_FILE = Resources.getResource("chord/v2/sample_chord_prediction.txt").getPath();
    private static final double EPSILON = 1.0E-10;

    private static final double V1_BRCA1 = 0.230;
    private static final double V1_BRCA2 = 0.400;
    private static final double V1_HRD = 0.630;

    private static final double V2_BRCA1 = 0.1;
    private static final double V2_BRCA2 = 0.2;
    private static final double V2_HRD = 0.3;
    private static final ChordStatus V2_HR_STATUS = ChordStatus.HR_PROFICIENT;
    private static final String V2_HRD_TYPE = "none";

    @Test
    public void canLoadV1ChordFile() throws IOException {
        ChordAnalysis chordAnalysis = ChordFileReader.read(V1_CHORD_FILE);

        assertEquals(V1_BRCA1, chordAnalysis.BRCA1Value(), EPSILON);
        assertEquals(V1_BRCA2, chordAnalysis.BRCA2Value(), EPSILON);
        assertEquals(V1_HRD, chordAnalysis.hrdValue(), EPSILON);
        assertEquals(ChordStatus.UNKNOWN, chordAnalysis.hrStatus());
        assertEquals(ChordFileReader.V1_NA, chordAnalysis.hrdType());
        assertEquals(ChordFileReader.V1_NA, chordAnalysis.remarksHrStatus());
        assertEquals(ChordFileReader.V1_NA, chordAnalysis.remarksHrdType());
    }

    @Test
    public void canLoadV2ChordFile() throws IOException {
        ChordAnalysis chordAnalysis = ChordFileReader.read(V2_CHORD_FILE);

        assertEquals(V2_BRCA1, chordAnalysis.BRCA1Value(), EPSILON);
        assertEquals(V2_BRCA2, chordAnalysis.BRCA2Value(), EPSILON);
        assertEquals(V2_HRD, chordAnalysis.hrdValue(), EPSILON);
        assertEquals(V2_HR_STATUS, chordAnalysis.hrStatus());
        assertEquals(V2_HRD_TYPE, chordAnalysis.hrdType());
        assertEquals(Strings.EMPTY, chordAnalysis.remarksHrStatus());
        assertEquals(Strings.EMPTY, chordAnalysis.remarksHrdType());
    }

    @Test
    public void canConvertHRDTOStatus() {
        assertEquals(ChordStatus.CANNOT_BE_DETERMINED, ChordFileReader.extractHrStatus("cannot_be_determined"));
        assertEquals(ChordStatus.HR_PROFICIENT, ChordFileReader.extractHrStatus("HR_proficient"));
        assertEquals(ChordStatus.HR_DEFICIENT, ChordFileReader.extractHrStatus("HR_deficient"));
        assertEquals(ChordStatus.UNKNOWN, ChordFileReader.extractHrStatus("dgdfg"));
    }
}