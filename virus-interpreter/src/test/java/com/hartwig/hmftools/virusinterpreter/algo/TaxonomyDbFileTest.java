package com.hartwig.hmftools.virusinterpreter.algo;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class TaxonomyDbFileTest {

    private static final String TAXONOMY_DB_TSV = Resources.getResource("virus_interpreter/taxonomy_db.tsv").getPath();

    @Test
    public void canReadTaxonomyDbTsv() throws IOException {
        TaxonomyDb taxonomyDb = TaxonomyDbFile.loadFromTsv(TAXONOMY_DB_TSV);
        assertEquals(1, taxonomyDb.count());
    }
}