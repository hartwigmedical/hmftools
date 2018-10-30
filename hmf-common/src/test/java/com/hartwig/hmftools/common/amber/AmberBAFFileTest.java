package com.hartwig.hmftools.common.amber;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.junit.Test;

public class AmberBAFFileTest {

    private static final String OLD_AMBER_BAF_PATH = Resources.getResource("amber/old.amber.baf").getPath();
    private static final String NEW_AMBER_BAF_PATH = Resources.getResource("amber/new.amber.baf").getPath();

    private static final double EPSILON = 1e-4;

    @Test
    public void testNewOldCompatibility() throws IOException {
        final List<AmberBAF> oldFormat = Lists.newArrayList(AmberBAFFile.read(OLD_AMBER_BAF_PATH).get(HumanChromosome._1));
        final List<AmberBAF> newFormat = Lists.newArrayList(AmberBAFFile.read(NEW_AMBER_BAF_PATH).get(HumanChromosome._1));

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
