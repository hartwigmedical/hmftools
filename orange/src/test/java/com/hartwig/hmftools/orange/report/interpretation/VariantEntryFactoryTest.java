package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

import org.junit.Test;

public class VariantEntryFactoryTest {

    @Test
    public void canGenerateVariantEvents() {
        assertEquals("p.Gly12Cys",
                VariantEntryFactory.toVariantEvent("p.Gly12Cys",
                        "c.123A>C",
                        Sets.newHashSet(VariantEffect.MISSENSE),
                        CodingEffect.MISSENSE));
        assertEquals("c.123A>C splice",
                VariantEntryFactory.toVariantEvent("p.?", "c.123A>C", Sets.newHashSet(VariantEffect.MISSENSE), CodingEffect.SPLICE));
        assertEquals("c.123A>C",
                VariantEntryFactory.toVariantEvent("", "c.123A>C", Sets.newHashSet(VariantEffect.MISSENSE), CodingEffect.MISSENSE));
        assertEquals("missense_variant",
                VariantEntryFactory.toVariantEvent("", "", Sets.newHashSet(VariantEffect.MISSENSE), CodingEffect.MISSENSE));
    }
}