package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.*;

import org.junit.Test;

public class LimsGermlineFindingsChoiceTest {

    @Test
    public void canExtractGermlineChoice() {
        assertEquals(LimsGermlineFindingsChoice.ALL_ACTIONABLE,
                LimsGermlineFindingsChoice.extractChoiceInformedConsent("1: Behandelbare toevalsbevindingen", "WIDE02991111T"));
        assertEquals(LimsGermlineFindingsChoice.ALL,
                LimsGermlineFindingsChoice.extractChoiceInformedConsent("2: Alle toevalsbevindingen", "WIDE02991111T"));
        assertEquals(LimsGermlineFindingsChoice.NONE_FAMILY,
                LimsGermlineFindingsChoice.extractChoiceInformedConsent("3: Geen toevalsbevindingen, familie mag deze wel opvragen",
                        "WIDE02991111T"));
        assertEquals(LimsGermlineFindingsChoice.NONE,
                LimsGermlineFindingsChoice.extractChoiceInformedConsent("4: Geen toevalsbevindingen, familie mag deze niet opvragen",
                        "WIDE02991111T"));
        assertEquals(LimsGermlineFindingsChoice.UNKNOWN, LimsGermlineFindingsChoice.extractChoiceInformedConsent("", "CPCT02991111T"));
        assertEquals(LimsGermlineFindingsChoice.UNKNOWN, LimsGermlineFindingsChoice.extractChoiceInformedConsent("", "DRUP02991111T"));
        assertEquals(LimsGermlineFindingsChoice.UNKNOWN, LimsGermlineFindingsChoice.extractChoiceInformedConsent("", "COLO02991111T"));
        assertEquals(LimsGermlineFindingsChoice.UNKNOWN, LimsGermlineFindingsChoice.extractChoiceInformedConsent("", "CORE02991111T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoice() {
        LimsGermlineFindingsChoice.extractChoiceInformedConsent("ALL", "WIDE02991111T");
    }

}