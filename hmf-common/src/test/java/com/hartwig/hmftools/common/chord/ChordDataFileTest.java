package com.hartwig.hmftools.common.chord;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ChordDataFileTest
{
    private static final String CHORD_FILE = Resources.getResource("chord/sample_chord_prediction.txt").getPath();
    private static final double EPSILON = 1.0E-10;

    private static final double BRCA1 = 0.1;
    private static final double BRCA2 = 0.2;
    private static final double HRD = 0.3;
    private static final ChordStatus HR_STATUS = ChordStatus.HR_PROFICIENT;
    private static final String HRD_TYPE = "none";

    @Test
    public void canLoadChordFile() throws IOException
    {
        ChordData chordData = ChordDataFile.read(CHORD_FILE);

        assertEquals(BRCA1, chordData.BRCA1Value(), EPSILON);
        assertEquals(BRCA2, chordData.BRCA2Value(), EPSILON);
        assertEquals(HRD, chordData.hrdValue(), EPSILON);
        assertEquals(HR_STATUS, chordData.hrStatus());
        assertEquals(HRD_TYPE, chordData.hrdType());
        assertEquals(Strings.EMPTY, chordData.remarksHrStatus());
        assertEquals(Strings.EMPTY, chordData.remarksHrdType());
    }

    @Test
    public void canConvertHRDToStatus()
    {
        assertEquals(ChordStatus.CANNOT_BE_DETERMINED, ChordDataFile.extractHrStatus("cannot_be_determined"));
        assertEquals(ChordStatus.HR_PROFICIENT, ChordDataFile.extractHrStatus("HR_proficient"));
        assertEquals(ChordStatus.HR_DEFICIENT, ChordDataFile.extractHrStatus("HR_deficient"));
        assertEquals(ChordStatus.UNKNOWN, ChordDataFile.extractHrStatus("dgdfg"));
    }
}