package com.hartwig.hmftools.serve.transvar;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class TransvarCuratorTest {

    @Test
    public void canCurateExceptionalGenes() {
        TransvarCurator curator = new TransvarCurator();

        assertEquals("C17orf96", curator.curateGene("EPOP"));
        assertEquals("CCDC186", curator.curateGene("CCDC186"));
    }

    @Test
    public void canCurateProteinAnnotation() {
        TransvarCurator curator = new TransvarCurator();
        assertEquals("M1I", curator.curateProteinAnnotation("M1?"));
        assertEquals("V600E", curator.curateGene("V600E"));
    }
}