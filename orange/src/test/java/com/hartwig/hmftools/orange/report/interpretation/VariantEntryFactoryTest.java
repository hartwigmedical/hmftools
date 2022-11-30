package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VariantEntryFactoryTest {

    @Test
    public void canDetermineTranscriptImpact() {
        assertEquals("p.Gly12Cys",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact("c.123A>C")
                        .hgvsProteinImpact("p.Gly12Cys")
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.MISSENSE)
                        .build()));
        assertEquals("c.123A>C splice",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact("c.123A>C")
                        .hgvsProteinImpact("p.?")
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.SPLICE)
                        .build()));
        assertEquals("c.123A>C",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact("c.123A>C")
                        .hgvsProteinImpact(Strings.EMPTY)
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.MISSENSE)
                        .build()));
        assertEquals("missense_variant",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact(Strings.EMPTY)
                        .hgvsProteinImpact(Strings.EMPTY)
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.MISSENSE)
                        .build()));
    }
}