package com.hartwig.hmftools.orange.cohort.mapping;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CohortMappingFileTest {

    private static final String EXAMPLE_MAPPING_TSV = Resources.getResource("cohort/mapping/example_cohort_mapping.tsv").getPath();

    @Test
    public void canReadExampleCohortMappingTsv() throws IOException {
        List<CohortMapping> mappings = CohortMappingFile.read(EXAMPLE_MAPPING_TSV);
        assertEquals(7, mappings.size());

        CohortMapping lymphoid = find(mappings, "Lymphoid tissue");
        assertEquals(1, lymphoid.include().size());
        assertEquals(0, lymphoid.exclude().size());

        CohortMapping ovary = find(mappings, "Ovary");
        assertEquals(2, ovary.include().size());
        assertEquals(9, ovary.exclude().size());
    }

    @NotNull
    private static CohortMapping find(@NotNull List<CohortMapping> mappings, @NotNull String cancerType) {
        for (CohortMapping mapping : mappings) {
            if (mapping.cancerType().equals(cancerType)) {
                return mapping;
            }
        }

        throw new IllegalStateException("Could not find mapping for cancer type: " + cancerType);
    }

}