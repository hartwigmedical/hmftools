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
    public void testNullOrEmptyInputReturnsNull()
    {
        final LiftBackResolver resolver = newResolver();
        assertNull(SaTagRewriter.rewriteSaTag(null, resolver));
        assertNull(SaTagRewriter.rewriteSaTag("", resolver));
    }

    @Test
    public void testRefContigPassesThrough()
    {
        final LiftBackResolver resolver = newResolver();
        final String sa =CHR_1 + ",1000,+,50M,60,2;";
        final String rewritten =SaTagRewriter.rewriteSaTag(sa, resolver);
        assertEquals(sa, rewritten);
    }

    @Test
    public void testTxContigLifted()
    {
        final LiftBackResolver resolver = newResolver();
        final String sa =TX_CONTIG + ",1,+,50M,60,2;";
        final String rewritten =SaTagRewriter.rewriteSaTag(sa, resolver);
        assertEquals(CHR_1 + ",100,+,50M,60,2;", rewritten);
    }

    @Test
    public void testMalformedEntriesSkipped()
    {
        final LiftBackResolver resolver = newResolver();
        final String sa ="junk;"
                + CHR_1 + ",notanumber,+,50M,60,0;"
                + CHR_1 + ",1000,+,50M,60,2;";
        final String rewritten =SaTagRewriter.rewriteSaTag(sa, resolver);
        assertEquals(CHR_1 + ",1000,+,50M,60,2;", rewritten);
    }

    @Test
    public void testDuplicateLiftedEntriesDeduped()
    {
        final LiftBackResolver resolver = newResolver();
        final String entry = CHR_1 + ",1000,+,50M,60,2;";
        final String rewritten =SaTagRewriter.rewriteSaTag(entry + entry, resolver);
        assertEquals(entry, rewritten);
    }

    @Test
    public void testAllEntriesFailReturnsNull()
    {
        final LiftBackResolver resolver = newResolver();
        final String sa =TX_CONTIG + ",10000,+,50M,60,0;malformed;";
        assertNull(SaTagRewriter.rewriteSaTag(sa, resolver));
    }
}
