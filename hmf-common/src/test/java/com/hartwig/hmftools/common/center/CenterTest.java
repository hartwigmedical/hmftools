package com.hartwig.hmftools.common.center;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class CenterTest {

    private static final String CENTER_RESOURCE = Resources.getResource("center/centers.csv").getPath();

    @Test
    public void canReadFromFile() throws IOException {
        final CenterModel centerModel = Center.readFromCSV(CENTER_RESOURCE);
        assertEquals(2, centerModel.centerPerId().size());
        assertEquals(2, centerModel.centerPerHospital().size());
    }
}
