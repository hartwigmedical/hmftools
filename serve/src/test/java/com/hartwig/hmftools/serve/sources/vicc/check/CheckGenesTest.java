package com.hartwig.hmftools.serve.sources.vicc.check;

import static org.junit.Assert.*;

import org.junit.Test;

public class CheckGenesTest {

    @Test
    public void canDetermineUsingGenesForCuration() {
        assertTrue(CheckGenes.checkGensInPanelForCuration("IGK", "IGK-MYC Fusion"));
        assertFalse(CheckGenes.checkGensInPanelForCuration("BRAF", "PAPSS1-BRAF"));

    }

}