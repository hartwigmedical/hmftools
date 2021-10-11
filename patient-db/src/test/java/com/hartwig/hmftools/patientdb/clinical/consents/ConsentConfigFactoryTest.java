package com.hartwig.hmftools.patientdb.clinical.consents;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class ConsentConfigFactoryTest {

    private static final String INFORMED_CONSENTS_TSV = Resources.getResource("consents/informed_consents.tsv").getPath();

    @Test
    public void canReadConsentConfigFile() throws IOException {
        ConsentConfig consentConfig1 = ConsentConfigFactory.read(INFORMED_CONSENTS_TSV).get("1");

        assertEquals("1", consentConfig1.pifVersion());
        assertEquals(Lists.newArrayList("CPCT", "CPCT1"), consentConfig1.cohort());
        assertNull(consentConfig1.inHMF());
        assertNull(consentConfig1.outsideEU());
        assertEquals("",consentConfig1.pif222());
        assertNull(consentConfig1.pif222Values());
        assertEquals("",consentConfig1.pif221());
        assertNull(consentConfig1.pif221Values());
        assertEquals("Ja",consentConfig1.pif26HMF());
        assertEquals(Lists.newArrayList("Ja", "Nee"), consentConfig1.pif26HMFValues());
        assertEquals("Ja",consentConfig1.pif26BUG());
        assertEquals(Lists.newArrayList("Ja", "Nee"), consentConfig1.pif26BUGValues());

        ConsentConfig consentConfig2 = ConsentConfigFactory.read(INFORMED_CONSENTS_TSV).get("2");
        assertEquals("2", consentConfig2.pifVersion());
        assertEquals(Lists.newArrayList("CPCT"), consentConfig2.cohort());
        assertNull(consentConfig2.inHMF());
        assertNull(consentConfig2.outsideEU());
        assertEquals("Yes",consentConfig2.pif222());
        assertEquals(Lists.newArrayList("Yes", "No"), consentConfig2.pif222Values());
        assertEquals("Yes",consentConfig2.pif221());
        assertEquals(Lists.newArrayList("Yes", "No"), consentConfig2.pif221Values());
        assertNull(consentConfig2.pif26HMF());
        assertNull(consentConfig2.pif26HMFValues());
        assertNull(consentConfig2.pif26BUG());
        assertNull(consentConfig2.pif26BUGValues());

        ConsentConfig consentConfig3 = ConsentConfigFactory.read(INFORMED_CONSENTS_TSV).get("3");
        assertEquals("3", consentConfig3.pifVersion());
        assertEquals(Lists.newArrayList("WIDE"), consentConfig3.cohort());
        assertTrue(consentConfig3.inHMF());
        assertTrue(consentConfig3.outsideEU());
        assertNull(consentConfig3.pif222());
        assertNull(consentConfig3.pif222Values());
        assertNull(consentConfig3.pif221());
        assertNull(consentConfig3.pif221Values());
        assertNull(consentConfig3.pif26HMF());
        assertNull(consentConfig3.pif26HMFValues());
        assertNull(consentConfig3.pif26BUG());
        assertNull(consentConfig3.pif26BUGValues());

        ConsentConfig consentConfig4 = ConsentConfigFactory.read(INFORMED_CONSENTS_TSV).get("4");
        assertEquals("4", consentConfig4.pifVersion());
        assertEquals(Lists.newArrayList("CPCT"), consentConfig4.cohort());
        assertNull(consentConfig4.inHMF());
        assertNull(consentConfig4.outsideEU());
        assertNull(consentConfig4.pif222());
        assertNull(consentConfig4.pif222Values());
        assertNull(consentConfig4.pif221());
        assertNull(consentConfig4.pif221Values());
        assertNull(consentConfig4.pif26HMF());
        assertNull(consentConfig4.pif26HMFValues());
        assertNull(consentConfig4.pif26BUG());
        assertNull(consentConfig4.pif26BUGValues());
    }
}