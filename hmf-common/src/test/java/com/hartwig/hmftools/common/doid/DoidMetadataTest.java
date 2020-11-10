package com.hartwig.hmftools.common.doid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class DoidMetadataTest {

    @Test
    public void canResolveSnomedID() {
        DoidMetadata emptyMetadata = ImmutableDoidMetadata.builder().build();
        assertNull(emptyMetadata.snomedId());

        DoidMetadata properMetadata = ImmutableDoidMetadata.builder()
                .addXrefs(ImmutableDoidXref.builder().val("DoidXref{val=SNOMEDCT_US_2020_03_01:109355002}").build())
                .build();

        assertEquals("109355002", properMetadata.snomedId());

        DoidMetadata differentID = ImmutableDoidMetadata.builder()
                .addXrefs(ImmutableDoidXref.builder().val("DoidXref{val=HMF:1}").build())
                .build();

        assertNull(differentID.snomedId());

        DoidMetadata wrongSnomedFormat = ImmutableDoidMetadata.builder()
                .addXrefs(ImmutableDoidXref.builder().val("SNOMED found!").build())
                .build();

        assertNull(wrongSnomedFormat.snomedId());
    }
}