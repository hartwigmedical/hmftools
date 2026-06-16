package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import org.junit.Test;

public class SaTagRewriterTest
{
    private static LiftBackResolver newResolver()
    {
        return new LiftBackResolver(List.of(threeExonContig()));
    }

    @Test
    public void testRewriteSaTag()
    {
        final LiftBackResolver resolver = newResolver();

        // null or empty input -> null
        assertNull(SaTagRewriter.rewriteSaTag(null, resolver));
        assertNull(SaTagRewriter.rewriteSaTag("", resolver));

        // ref-contig entry passes through unchanged
        final String refEntry = CHR_1 + ",1000,+,50M,60,2;";
        assertEquals(refEntry, SaTagRewriter.rewriteSaTag(refEntry, resolver));

        // tx-contig entry lifted to genomic coordinates
        assertEquals(CHR_1 + ",100,+,50M,60,2;", SaTagRewriter.rewriteSaTag(TX_CONTIG + ",1,+,50M,60,2;", resolver));

        // malformed entries skipped, valid one survives
        assertEquals(refEntry, SaTagRewriter.rewriteSaTag(
                "junk;" + CHR_1 + ",notanumber,+,50M,60,0;" + refEntry, resolver));

        // duplicate lifted entries deduped
        assertEquals(refEntry, SaTagRewriter.rewriteSaTag(refEntry + refEntry, resolver));

        // all entries fail -> null
        assertNull(SaTagRewriter.rewriteSaTag(TX_CONTIG + ",10000,+,50M,60,0;malformed;", resolver));
    }
}
