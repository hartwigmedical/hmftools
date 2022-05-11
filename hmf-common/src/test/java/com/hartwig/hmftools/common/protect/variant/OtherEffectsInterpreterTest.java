package com.hartwig.hmftools.common.protect.variant;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.junit.Test;

public class OtherEffectsInterpreterTest {

    @Test
    public void canExtractFromOtherEffects() {
        String example = "ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT";

        assertEquals("ENST00000579755", OtherEffectsInterpreter.transcript(example));
        assertEquals("c.246_247delCG", OtherEffectsInterpreter.hgvsCodingImpact(example));
        assertEquals("p.Gly83fs", OtherEffectsInterpreter.hgvsProteinImpact(example));
        assertEquals("frameshift_variant", OtherEffectsInterpreter.effect(example));
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, OtherEffectsInterpreter.codingEffect(example));
    }
}