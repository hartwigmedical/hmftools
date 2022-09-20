package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;

import org.apache.logging.log4j.util.Strings;
import org.junit.Assert;
import org.junit.Test;

public class AltTranscriptReportableInfoTest
{
    @Test
    public void canExtractFromOtherEffects() {
        String example = "ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT";

        Assert.assertEquals("ENST00000579755", AltTranscriptReportableInfo.firstOtherTranscript(example));
        assertEquals("c.246_247delCG", AltTranscriptReportableInfo.firstOtherHgvsCodingImpact(example));
        assertEquals("p.Gly83fs", AltTranscriptReportableInfo.firstOtherHgvsProteinImpact(example));
        assertEquals("frameshift_variant", AltTranscriptReportableInfo.firstOtherEffects(example));
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, AltTranscriptReportableInfo.firstOtherCodingEffect(example));
    }

    @Test
    public void worksOnEmptyString() {
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherTranscript(Strings.EMPTY));
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherHgvsCodingImpact(Strings.EMPTY));
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherHgvsProteinImpact(Strings.EMPTY));
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherEffects(Strings.EMPTY));
        assertEquals(CodingEffect.UNDEFINED, AltTranscriptReportableInfo.firstOtherCodingEffect(Strings.EMPTY));
    }
}