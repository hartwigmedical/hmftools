package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.*;

import org.junit.Test;

public class LimsInformedConsentTest {

    @Test
    public void canExtractGermlineChoice() {
        assertEquals(LimsInformedConsent.ALL_ACTIONABLE,
                LimsInformedConsent.extractChoiceInformedConsent("1: Behandelbare toevalsbevindingen", "WIDE02991111T"));
        assertEquals(LimsInformedConsent.ALL,
                LimsInformedConsent.extractChoiceInformedConsent("2: Alle toevalsbevindingen", "WIDE02991111T"));
        assertEquals(LimsInformedConsent.NONE_FAMILY,
                LimsInformedConsent.extractChoiceInformedConsent("3: Geen toevalsbevindingen, familie mag deze wel opvragen",
                        "WIDE02991111T"));
        assertEquals(LimsInformedConsent.NONE,
                LimsInformedConsent.extractChoiceInformedConsent("4: Geen toevalsbevindingen, familie mag deze niet opvragen",
                        "WIDE02991111T"));
        assertEquals(LimsInformedConsent.UNKNOWN, LimsInformedConsent.extractChoiceInformedConsent("", "CPCT02991111T"));
        assertEquals(LimsInformedConsent.UNKNOWN, LimsInformedConsent.extractChoiceInformedConsent("", "DRUP02991111T"));
        assertEquals(LimsInformedConsent.UNKNOWN, LimsInformedConsent.extractChoiceInformedConsent("", "COLO02991111T"));
        assertEquals(LimsInformedConsent.UNKNOWN, LimsInformedConsent.extractChoiceInformedConsent("", "CORE02991111T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoice() {
        LimsInformedConsent.extractChoiceInformedConsent("ALL", "WIDE02991111T");
    }

}