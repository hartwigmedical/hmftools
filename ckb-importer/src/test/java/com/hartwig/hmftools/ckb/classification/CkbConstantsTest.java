package com.hartwig.hmftools.ckb.classification;

import static org.junit.Assert.assertFalse;

import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping37to38;

import org.junit.Test;

public class CkbConstantsTest {

    @Test
    public void unresolvableAndUnmappableGenesDoNotExistIn38() {
        GeneNameMapping37to38 geneNameMapping = GeneNameMapping37to38.loadFromEmbeddedResource();

        for (String unresolvableGene : CkbConstants.NON_EXISTING_GENES) {
            assertFalse(geneNameMapping.isValidV38Gene(unresolvableGene));
        }

        for (String unmappableGene : CkbConstants.UNMAPPABLE_GENES) {
            assertFalse(geneNameMapping.isValidV38Gene(unmappableGene));
        }
    }
}