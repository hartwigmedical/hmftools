package com.hartwig.hmftools.common.pcf;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PCFFileTest {

    private static final int WINDOW = 1000;
    private static final String BASE_PATH = Resources.getResource("pcf").getPath() + File.separator;

    @Test
    public void testBafFile() throws IOException {

        final ListMultimap<String, PCFPosition> resultMap =
                PCFFile.readPositions(WINDOW, PCFSource.TUMOR_BAF, BASE_PATH + File.separator + "baf.pcf");

        final List<PCFPosition> chromOneResults = resultMap.get("1");

        assertEquals(4, chromOneResults.size());
        assertPosition(93548001, chromOneResults.get(0));
        assertPosition(193800001, chromOneResults.get(1));
        assertPosition(193803001, chromOneResults.get(2));
        assertPosition(193804001, chromOneResults.get(3));

        final List<PCFPosition> chromTwoResults = resultMap.get("2");

        assertEquals(4, chromTwoResults.size());
        assertPosition(93548001, chromTwoResults.get(0));
        assertPosition(193800001, chromTwoResults.get(1));
        assertPosition(193803001, chromTwoResults.get(2));
        assertPosition(193804001, chromTwoResults.get(3));
    }

    @Test
    public void testRatioFile() throws IOException {
        final List<PCFPosition> results =
                PCFFile.readPositions(WINDOW, PCFSource.TUMOR_BAF, BASE_PATH + File.separator + "ratio.pcf").get("1");

        assertEquals(5, results.size());
        assertPosition(835001, results.get(0));
        assertPosition(2583001, results.get(1));
        assertPosition(2584001, results.get(2));
        assertPosition(2695001, results.get(3));
        assertPosition(4363001, results.get(4));
    }

    private static void assertPosition(long position, @NotNull final PCFPosition victim) {
        assertEquals(position, victim.position());
    }
}
