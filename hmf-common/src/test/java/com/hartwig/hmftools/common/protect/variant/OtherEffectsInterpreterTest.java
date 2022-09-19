package com.hartwig.hmftools.common.protect.variant;

import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.apache.logging.log4j.util.Strings;
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

    @Test
    public void worksWithSplicing() {
        String example = "ENST00000579755|c.194-3_194-2delCA|p.?|splice_acceptor_variant&intron_variant|SPLICE";
        assertEquals(CodingEffect.SPLICE, OtherEffectsInterpreter.codingEffect(example));
    }

    @Test
    public void worksOnEmptyString() {
        assertEquals(Strings.EMPTY, OtherEffectsInterpreter.transcript(Strings.EMPTY));
        assertEquals(Strings.EMPTY, OtherEffectsInterpreter.hgvsCodingImpact(Strings.EMPTY));
        assertEquals(Strings.EMPTY, OtherEffectsInterpreter.hgvsProteinImpact(Strings.EMPTY));
        assertEquals(Strings.EMPTY, OtherEffectsInterpreter.effect(Strings.EMPTY));
        assertEquals(CodingEffect.UNDEFINED, OtherEffectsInterpreter.codingEffect(Strings.EMPTY));
    }

    @Test
    public void worksOnDoubleExample() {
        String impact = "ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT";
        String example = impact+ VAR_IMPACT_OTHER_REPORT_DELIM + impact;

        assertEquals("ENST00000579755", OtherEffectsInterpreter.transcript(example));
        assertEquals("c.246_247delCG", OtherEffectsInterpreter.hgvsCodingImpact(example));
        assertEquals("p.Gly83fs", OtherEffectsInterpreter.hgvsProteinImpact(example));
        assertEquals("frameshift_variant", OtherEffectsInterpreter.effect(example));
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, OtherEffectsInterpreter.codingEffect(example));
    }
}