package com.hartwig.hmftools.common.chord;

import static org.junit.Assert.*;

import org.junit.Test;

public class ChordStatusTest {

    @Test
    public void canConvertHRDTOStatus() {
        assertEquals(ChordStatus.fromHRD(0.4), ChordStatus.HRP);
        assertEquals(ChordStatus.fromHRD(0.7), ChordStatus.HRD);
        assertEquals(ChordStatus.fromHRD(0.5), ChordStatus.HRD);
    }

}