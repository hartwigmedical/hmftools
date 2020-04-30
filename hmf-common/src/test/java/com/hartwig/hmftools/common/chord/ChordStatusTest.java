package com.hartwig.hmftools.common.chord;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ChordStatusTest {

    @Test
    public void canConvertHRDTOStatus() {
        assertEquals(ChordStatus.HRP, ChordStatus.fromHRD(0.4));
        assertEquals(ChordStatus.HRD, ChordStatus.fromHRD(0.7));
        assertEquals(ChordStatus.HRD, ChordStatus.fromHRD(0.5));
    }
}