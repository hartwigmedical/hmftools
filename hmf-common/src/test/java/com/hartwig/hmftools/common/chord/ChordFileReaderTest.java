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
    private static final String V2_HR_STATUS = "HR_proficient";
    private static final String V2_HRD_TYPE = "none";

    @Test
    public void canLoadV1ChordFile() throws IOException {
        ChordAnalysis chordAnalysis = ChordFileReader.read(V1_CHORD_FILE);

        assertEquals(V1_BRCA1, chordAnalysis.BRCA1Value(), EPSILON);
        assertEquals(V1_BRCA2, chordAnalysis.BRCA2Value(), EPSILON);
        assertEquals(V1_HRD, chordAnalysis.hrdValue(), EPSILON);
        assertEquals(ChordFileReader.extractHrStatus(ChordFileReader.V1_NA), chordAnalysis.hrStatus());
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
        assertEquals(ChordFileReader.extractHrStatus(V2_HR_STATUS), chordAnalysis.hrStatus());
        assertEquals(V2_HRD_TYPE, chordAnalysis.hrdType());
        assertEquals(Strings.EMPTY, chordAnalysis.remarksHrStatus());
        assertEquals(Strings.EMPTY, chordAnalysis.remarksHrdType());
    }

    @Test
    public void canConvertHRDTOStatus() {
        ChordAnalysis chordAnalysisCannotBeDetermined = ImmutableChordAnalysis.builder()
                .BRCA1Value(0.10)
                .BRCA2Value(0.20)
                .hrdValue(0.30)
                .hrStatus(ChordStatus.CANNOT_BE_DETERMINED)
                .hrdType(Strings.EMPTY)
                .remarksHrStatus(Strings.EMPTY)
                .remarksHrdType(Strings.EMPTY)
                .build();

        assertEquals(ChordFileReader.extractHrStatus("cannot_be_determined"), chordAnalysisCannotBeDetermined.hrStatus());

        ChordAnalysis chordAnalysisHRP = ImmutableChordAnalysis.builder()
                .BRCA1Value(0.10)
                .BRCA2Value(0.20)
                .hrdValue(0.30)
                .hrStatus(ChordStatus.HR_PROFICIENT)
                .hrdType(Strings.EMPTY)
                .remarksHrStatus(Strings.EMPTY)
                .remarksHrdType(Strings.EMPTY)
                .build();

        assertEquals(ChordFileReader.extractHrStatus("HR_proficient"), chordAnalysisHRP.hrStatus());

        ChordAnalysis chordAnalysisHRD = ImmutableChordAnalysis.builder()
                .BRCA1Value(0.10)
                .BRCA2Value(0.20)
                .hrdValue(0.30)
                .hrStatus(ChordStatus.HR_DEFICIENT)
                .hrdType(Strings.EMPTY)
                .remarksHrStatus(Strings.EMPTY)
                .remarksHrdType(Strings.EMPTY)
                .build();

        assertEquals(ChordFileReader.extractHrStatus("HR_deficient"), chordAnalysisHRD.hrStatus());

        ChordAnalysis chordAnalysisUnknown = ImmutableChordAnalysis.builder()
                .BRCA1Value(0.10)
                .BRCA2Value(0.20)
                .hrdValue(0.30)
                .hrStatus(ChordStatus.UNKNOWN)
                .hrdType(Strings.EMPTY)
                .remarksHrStatus(Strings.EMPTY)
                .remarksHrdType(Strings.EMPTY)
                .build();

        assertEquals(ChordFileReader.extractHrStatus("dgdfg"), chordAnalysisUnknown.hrStatus());

    }
}