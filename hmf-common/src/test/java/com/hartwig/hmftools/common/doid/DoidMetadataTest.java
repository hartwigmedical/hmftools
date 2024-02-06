package com.hartwig.hmftools.common.doid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class DoidMetadataTest
{
    @Test
    public void canResolveSnomedConceptID()
    {
        DoidMetadata emptyMetadata = ImmutableDoidMetadata.builder().build();
        assertNull(emptyMetadata.snomedConceptId());

        DoidMetadata properMetadata = ImmutableDoidMetadata.builder()
                .addXrefs(ImmutableDoidXref.builder().val("SNOMEDCT_US_2020_03_01:109355002").build())
                .build();

        assertEquals("109355002", properMetadata.snomedConceptId());

        DoidMetadata differentID = ImmutableDoidMetadata.builder()
                .addXrefs(ImmutableDoidXref.builder().val("HMF:1").build())
                .build();

        assertNull(differentID.snomedConceptId());

        DoidMetadata wrongSnomedFormat = ImmutableDoidMetadata.builder()
                .addXrefs(ImmutableDoidXref.builder().val("SNOMED found!").build())
                .build();

        assertNull(wrongSnomedFormat.snomedConceptId());
    }
}