package com.hartwig.hmftools.serve.sources.vicc.check;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.junit.Ignore;
import org.junit.Test;

public class GeneCheckerTest {

    @Test
    public void canCorrectlyCheckGenes() {
        GeneChecker geneChecker = new GeneChecker();
        HmfTranscriptRegion canonicalTranscript = HmfGenePanelSupplier.allGenesMap37().get("BRAF");
        assertTrue(geneChecker.isValidGene("BRAF", canonicalTranscript, "Amplification", null));

        assertFalse(geneChecker.isValidGene("IGH", null, "Deletion", null));
        assertTrue(geneChecker.isValidGene("IGH", null, "IGH-MYC", "fusion"));

        assertFalse(geneChecker.isValidGene("I am not a gene", null, "event", null));

    }
}