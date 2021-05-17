package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class TaxonomyDbFileTest {

    private static final String TAXONOMY_DB_TSV = Resources.getResource("viral_reporting/taxonomy_db.tsv").getPath();

    @Test
    public void canReadTaxonomyDbTsv() throws IOException {
        TaxonomyDb taxonomyDb = TaxonomyDbFile.loadFromTsv(TAXONOMY_DB_TSV);
        assertEquals(1, taxonomyDb.count());

        assertTrue(taxonomyDb.taxidExists(1));
        assertFalse(taxonomyDb.taxidExists(2));
    }
}