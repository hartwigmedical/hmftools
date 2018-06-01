package com.hartwig.hmftools.patientreporter.civic;

import static org.junit.Assert.assertNotNull;

import org.junit.Test;

public class AdditionalCivicMatchesTest {

    @Test
    public void canCreateAdditionalVariantsMapping() {
        assertNotNull(AdditionalCivicMatches.createAdditionalVariantsMapping());
    }
}