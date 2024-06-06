package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.parseAltTranscriptInfo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;

import org.apache.logging.log4j.util.Strings;
import org.junit.Assert;
import org.junit.Test;

public class AltTranscriptReportableInfoTest
{
    @Test
    public void canExtractFromOtherEffects()
    {
        String example = "ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT";

        Assert.assertEquals("ENST00000579755", AltTranscriptReportableInfo.firstOtherTranscript(example));
        assertEquals("c.246_247delCG", AltTranscriptReportableInfo.firstOtherHgvsCodingImpact(example));
        assertEquals("p.Gly83fs", AltTranscriptReportableInfo.firstOtherHgvsProteinImpact(example));
        assertEquals("frameshift_variant", AltTranscriptReportableInfo.firstOtherEffects(example));
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, AltTranscriptReportableInfo.firstOtherCodingEffect(example));
    }

    @Test
    public void canExtractFromMultipleOtherEffects()
    {
        AltTranscriptReportableInfo altInfo1 = new AltTranscriptReportableInfo(
                "CDKN2A", "ENST00000579755", "c.194-3_194-2delCA", "p.?",
                "splice_acceptor_variant&intron_variant", CodingEffect.SPLICE);

        AltTranscriptReportableInfo altInfo2 = new AltTranscriptReportableInfo(
                "SOME_GENE","ENST00000123456", "c.194-3_194-2delCA", "p.?",
                "splice_acceptor_variant&intron_variant", CodingEffect.NONSENSE_OR_FRAMESHIFT);

        String example = altInfo1.serialise() + VAR_IMPACT_OTHER_REPORT_DELIM + altInfo2.serialise();

        List<AltTranscriptReportableInfo> altTransInfos = parseAltTranscriptInfo(example);
        assertEquals(2, altTransInfos.size());
        assertTrue(altTransInfos.get(0).TransName.equals(altInfo1.TransName));
        assertEquals(altTransInfos.get(0).Effect, altInfo1.Effect);

        assertTrue(altTransInfos.get(1).TransName.equals(altInfo2.TransName));
        assertEquals(altTransInfos.get(1).Effect, altInfo2.Effect);
    }

    @Test
    public void worksOnEmptyString()
    {
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherTranscript(Strings.EMPTY));
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherHgvsCodingImpact(Strings.EMPTY));
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherHgvsProteinImpact(Strings.EMPTY));
        assertEquals(Strings.EMPTY, AltTranscriptReportableInfo.firstOtherEffects(Strings.EMPTY));
        assertEquals(CodingEffect.UNDEFINED, AltTranscriptReportableInfo.firstOtherCodingEffect(Strings.EMPTY));
    }
}