package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.HG19;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RefGenomeFunctionsTest {

    @Test
    public void canConvertChromosomes() {
        String chr37 = "10";
        String chr19 = "chr10";
        String chr38 = "chr10";

        assertEquals(chr37, RefGenomeFunctions.refGenomeChromosome(chr37, V37));
        assertEquals(chr38, RefGenomeFunctions.refGenomeChromosome(chr37, V38));

        assertEquals(chr19, RefGenomeFunctions.refGenomeChromosome(chr37, HG19));
        assertEquals(chr38, RefGenomeFunctions.refGenomeChromosome(chr38, HG19));

        assertEquals(chr37, RefGenomeFunctions.refGenomeChromosome(chr19, V37));
        assertEquals(chr37, RefGenomeFunctions.refGenomeChromosome(chr38, V37));
    }
}