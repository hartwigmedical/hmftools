package com.hartwig.hmftools.serve.curation;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;

import org.junit.Test;

public class DoidLookupFactoryTest {

    private static final String EXAMPLE_TSV = Resources.getResource("curation/doid_mapping_example.tsv").getPath();

    @Test
    public void canCreateFromTestResource() throws IOException {
        DoidLookup doidLookup = DoidLookupFactory.buildFromMappingTsv(EXAMPLE_TSV);

        assertEquals(2, doidLookup.cancerTypeToDoidsMapping().size());
        assertEquals(Sets.newHashSet("0123", "0124"), doidLookup.cancerTypeToDoidsMapping().get("cancerA"));
        assertEquals(Sets.newHashSet("0125"), doidLookup.cancerTypeToDoidsMapping().get("cancerB"));
    }
}