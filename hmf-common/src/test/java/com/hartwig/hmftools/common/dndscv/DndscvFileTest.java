package com.hartwig.hmftools.common.dndscv;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class DndscvFileTest {

    private static final String BASE_PATH = Resources.getResource("dndscv").getPath() + File.separator;
    private static final double EPSILON = 1e-4;

    @Test
    public void testLongVersion() throws IOException {
        final List<Dndscv> dndscvs = DndscvFile.read(BASE_PATH + "dndscv_long.txt");
        assertEquals("TP53", dndscvs.get(0).gene());
        assertEquals(0.159325033617401, dndscvs.get(dndscvs.size() - 1).qScore(), EPSILON);
        assertEquals(0.000793056414223003, dndscvs.get(dndscvs.size() - 1).pScore(), EPSILON);
    }

    @Test
    public void testShortVersion() throws IOException {
        final List<Dndscv> dndscvs = DndscvFile.read(BASE_PATH + "dndscv_short.txt");
        assertEquals("KRAS", dndscvs.get(0).gene());
        assertEquals(0.00604243162028628, dndscvs.get(dndscvs.size() - 1).qScore(), EPSILON);
        assertEquals(2.97760443209727e-05, dndscvs.get(dndscvs.size() - 1).pScore(), EPSILON);
    }
}
