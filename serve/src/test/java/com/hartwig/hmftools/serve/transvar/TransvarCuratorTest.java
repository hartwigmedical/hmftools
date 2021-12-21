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
}