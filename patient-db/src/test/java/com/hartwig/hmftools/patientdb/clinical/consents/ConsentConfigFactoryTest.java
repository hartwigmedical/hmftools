package com.hartwig.hmftools.patientdb.clinical.consents;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class ConsentConfigFactoryTest {

    private static final String INFORMED_CONSENTS_TSV = Resources.getResource("consents/informed_consents.tsv").getPath();

    @Test
    public void canReadConsentConfigFile() throws IOException {
        ConsentConfig consentConfig = ConsentConfigFactory.read(INFORMED_CONSENTS_TSV).get("1");

        assertEquals("1", consentConfig.pifVersion());
        assertEquals("",consentConfig.pif222());
        assertNull(consentConfig.pif222Values());
        assertEquals("",consentConfig.pif221());
        assertNull(consentConfig.pif221Values());
        assertEquals("Ja",consentConfig.pif26HMF());
        assertEquals(Lists.newArrayList("Ja", "Nee"), consentConfig.pif26HMFValues());
        assertEquals("Ja",consentConfig.pif26BUG());
        assertEquals(Lists.newArrayList("Ja", "Nee"), consentConfig.pif26BUGValues());
    }
}