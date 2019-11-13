package com.hartwig.hmftools.common.cobalt;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.junit.Test;

public class CobaltRatioFileTest {

    private static final String HG19_PATH = Resources.getResource("cobalt/hg19.cobalt.tsv").getPath();
    private static final String HG38_PATH = Resources.getResource("cobalt/hg38.cobalt.tsv").getPath();

    @Test
    public void testHG38() throws IOException {
        final List<CobaltRatio> hg38 = Lists.newArrayList(CobaltRatioFile.read(HG38_PATH).get(HumanChromosome._1));
        assertEquals(5, hg38.size());
    }

    @Test
    public void testHG19() throws IOException {
        final List<CobaltRatio> hg19 = Lists.newArrayList(CobaltRatioFile.read(HG19_PATH).get(HumanChromosome._1));
        assertEquals(4, hg19.size());
    }
}
