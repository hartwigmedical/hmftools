package com.hartwig.hmftools.common.pharmacogenetics;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PGXGenotypeFileTest {

    private static final String FILE = Resources.getResource("pharmacogenetics/sample_genotype.txt").getPath();

    private static final String GENE = "GENE";
    private static final String HAPLOTYPE = "*1_HET";
    private static final String FUNCTION = "Reduced Function";
    private static final String LINKED_DRUGS = "Capecitabine";
    private static final String URL_PRESCRIPTION_INFO = "link";
    private static final String PANEL_VERSION = "panel";
    private static final String REPO_VERSION = "pilot";

    @Test
    public void loadPgxGenotypeFile() throws IOException {
        List<PGXGenotype> pgxGenotypes = PGXGenotypeFile.read(FILE);

        assertEquals(1, pgxGenotypes.size());

        PGXGenotype genotype1 = pgxGenotypes.get(0);
        assertEquals(GENE, genotype1.gene());
        assertEquals(HAPLOTYPE, genotype1.haplotype());
        assertEquals(FUNCTION, genotype1.function());
        assertEquals(LINKED_DRUGS, genotype1.linkedDrugs());
        assertEquals(URL_PRESCRIPTION_INFO, genotype1.urlPrescriptionInfo());
        assertEquals(PANEL_VERSION, genotype1.panelVersion());
        assertEquals(REPO_VERSION, genotype1.repoVersion());
    }
}