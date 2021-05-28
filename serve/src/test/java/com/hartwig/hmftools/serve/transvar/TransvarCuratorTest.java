package com.hartwig.hmftools.serve.transvar;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Test;

public class TransvarCuratorTest {

    @Test
    public void canCurateExceptionalGenes() {
        TransvarCurator curator = new TransvarCurator(RefGenomeVersion.V38);

        assertEquals("C17orf96", curator.curateGene("EPOP"));
        assertEquals("CCDC186", curator.curateGene("CCDC186"));
    }
}