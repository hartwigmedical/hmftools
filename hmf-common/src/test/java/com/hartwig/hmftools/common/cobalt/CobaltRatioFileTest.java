package com.hartwig.hmftools.common.cobalt;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class CobaltRatioFileTest
{
    @SuppressWarnings("UnstableApiUsage")
    private static final String V37_PATH = Resources.getResource("cobalt/cobalt.37.tsv").getPath();
    @SuppressWarnings("UnstableApiUsage")
    private static final String V38_PATH = Resources.getResource("cobalt/cobalt.38.tsv").getPath();

    @Test
    public void testV38() throws IOException
    {
        final List<CobaltRatio> v38 = CobaltRatioFile.read(V38_PATH).get(HumanChromosome._1);
        assertEquals(5, v38.size());
    }

    @Test
    public void testV37() throws IOException
    {
        final List<CobaltRatio> v37 = CobaltRatioFile.read(V37_PATH).get(HumanChromosome._1);
        assertEquals(4, v37.size());
    }
}
