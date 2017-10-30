package com.hartwig.hmftools.common.amber;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class AmberBAFFileTest {

    private static final String BASE_PATH = Resources.getResource("amber").getPath() + File.separator;
    private static final double EPSILON = 1e-4;

    @Test
    public void testNewOldCompatibility() throws IOException {
        final List<AmberBAF> oldFormat = Lists.newArrayList(AmberBAFFile.read(BASE_PATH + "old.amber.baf").get("1"));
        final List<AmberBAF> newFormat = Lists.newArrayList(AmberBAFFile.read(BASE_PATH + "new.amber.baf").get("1"));

        for (int i = 0; i < oldFormat.size(); i++) {

            final AmberBAF oldBaf = oldFormat.get(i);
            final AmberBAF newBaf = newFormat.get(i);

            assertEquals(oldBaf.chromosome(), newBaf.chromosome());
            assertEquals(oldBaf.position(), newBaf.position());
            assertEquals(oldBaf.tumorBAF(), newBaf.tumorBAF(), EPSILON);
            assertEquals(oldBaf.tumorModifiedBAF(), newBaf.tumorModifiedBAF(), EPSILON);
            assertEquals(0, oldBaf.normalDepth());
            assertEquals(0, oldBaf.tumorDepth());
        }
    }
}
