package com.hartwig.hmftools.rose.actionability;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class TypeAlterationTest {

    @Test
    public void canExtractType() {
        assertEquals(TypeAlteration.ACTIVATING_MUTATION, TypeAlteration.toType("ACTIVATING_MUTATION"));
        assertEquals(TypeAlteration.AMPLIFICATION, TypeAlteration.toType("AMPLIFICATION"));
        assertEquals(TypeAlteration.EXTRACELLULAR_DOMAIN_MUTATION, TypeAlteration.toType("EXTRACELLULAR_DOMAIN_MUTATION"));
        assertEquals(TypeAlteration.FUSION, TypeAlteration.toType("FUSION"));
        assertEquals(TypeAlteration.INACTIVATION, TypeAlteration.toType("INACTIVATION"));
        assertEquals(TypeAlteration.INTERNAL_DELETION, TypeAlteration.toType("INTERNAL_DELETION"));
        assertEquals(TypeAlteration.KINASE_DOMAIN_DUPLICATION, TypeAlteration.toType("KINASE_DOMAIN_DUPLICATION"));
        assertEquals(TypeAlteration.LOSS, TypeAlteration.toType("LOSS"));
        assertEquals(TypeAlteration.POSITIVE, TypeAlteration.toType("POSITIVE"));
        assertEquals(TypeAlteration.RESISTANCE_MUTATION, TypeAlteration.toType("RESISTANCE_MUTATION"));
        assertEquals(TypeAlteration.PURITY, TypeAlteration.toType("PURITY"));
        assertEquals(TypeAlteration.PURITY_UNRELIABLE, TypeAlteration.toType("PURITY_UNRELIABLE"));
        assertEquals(TypeAlteration.NO_ONCOGENIC, TypeAlteration.toType("NO_ONCOGENIC"));
        assertEquals(TypeAlteration.NO_ACTIONABLE, TypeAlteration.toType("NO_ACTIONABLE"));
        assertEquals(TypeAlteration.FINDINGS, TypeAlteration.toType("FINDINGS"));
        assertEquals(TypeAlteration.GERMLINE, TypeAlteration.toType("GERMLINE"));
        assertEquals(TypeAlteration.CUPPA, TypeAlteration.toType("CUPPA"));
        assertEquals(TypeAlteration.CUPPA_INCONCLUSIVE, TypeAlteration.toType("CUPPA_INCONCLUSIVE"));
        assertEquals(TypeAlteration.NO_HRD_CAUSE, TypeAlteration.toType("NO_HRD_CAUSE"));
        assertEquals(TypeAlteration.NOT_BIALLELIC, TypeAlteration.toType("NOT_BIALLELIC"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownType() {
        TypeAlteration.toType("ABC");
    }

}