package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class TaxonomyDbTest {

    @Test
    public void canMatchTaxidToName() {
        Map<Integer, String> taxidToNameMap = Maps.newHashMap();
        taxidToNameMap.put(1, "species1");
        TaxonomyDb taxonomyDb = new TaxonomyDb(taxidToNameMap);

        assertEquals("species1", taxonomyDb.lookupName(1));
        assertNotEquals("species2", taxonomyDb.lookupName(1));

        assertTrue(taxonomyDb.taxidExists(1));
        assertFalse(taxonomyDb.taxidExists(3));
    }
}