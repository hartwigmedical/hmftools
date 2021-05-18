package com.hartwig.hmftools.common.peach;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PeachGenotypeFileTest {

    private static final String PEACH_GENOTYPES_FILE = Resources.getResource("peach/sample_genotype.txt").getPath();

    private static final String GENE = "GENE";
    private static final String HAPLOTYPE = "*1_HET";
    private static final String FUNCTION = "Reduced Function";
    private static final String LINKED_DRUGS = "Capecitabine";
    private static final String URL_PRESCRIPTION_INFO = "link";
    private static final String PANEL_VERSION = "panel";
    private static final String REPO_VERSION = "pilot";

    @Test
    public void loadPeachGenotypeFile() throws IOException {
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(PEACH_GENOTYPES_FILE);

        assertEquals(1, peachGenotypes.size());

        PeachGenotype peachGenotype = peachGenotypes.get(0);
        assertEquals(GENE, peachGenotype.gene());
        assertEquals(HAPLOTYPE, peachGenotype.haplotype());
        assertEquals(FUNCTION, peachGenotype.function());
        assertEquals(LINKED_DRUGS, peachGenotype.linkedDrugs());
        assertEquals(URL_PRESCRIPTION_INFO, peachGenotype.urlPrescriptionInfo());
        assertEquals(PANEL_VERSION, peachGenotype.panelVersion());
        assertEquals(REPO_VERSION, peachGenotype.repoVersion());
    }
}