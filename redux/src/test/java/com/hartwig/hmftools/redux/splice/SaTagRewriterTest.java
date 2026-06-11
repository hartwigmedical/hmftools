package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

public class SaTagRewriterTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    private static LiftBackResolver newResolver()
    {
        ContigEntry entry = new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
        return new LiftBackResolver(List.of(entry));
    }

    @Test
    public void testNullOrEmptyInputReturnsNull()
    {
        LiftBackResolver resolver = newResolver();
        assertNull(SaTagRewriter.rewriteSaTag(null, resolver));
        assertNull(SaTagRewriter.rewriteSaTag("", resolver));
    }

    @Test
    public void testRefContigPassesThrough()
    {
        LiftBackResolver resolver = newResolver();
        String sa = CHR_1 + ",1000,+,50M,60,2;";
        String rewritten = SaTagRewriter.rewriteSaTag(sa, resolver);
        assertEquals(sa, rewritten);
    }

    @Test
    public void testTxContigLifted()
    {
        LiftBackResolver resolver = newResolver();
        String sa = TX_CONTIG + ",1,+,50M,60,2;";
        String rewritten = SaTagRewriter.rewriteSaTag(sa, resolver);
        assertEquals(CHR_1 + ",100,+,50M,60,2;", rewritten);
    }

    @Test
    public void testMalformedEntriesSkipped()
    {
        LiftBackResolver resolver = newResolver();
        String sa = "junk;"
                + CHR_1 + ",notanumber,+,50M,60,0;"
                + CHR_1 + ",1000,+,50M,60,2;";
        String rewritten = SaTagRewriter.rewriteSaTag(sa, resolver);
        assertEquals(CHR_1 + ",1000,+,50M,60,2;", rewritten);
    }

    @Test
    public void testDuplicateLiftedEntriesDeduped()
    {
        LiftBackResolver resolver = newResolver();
        String entry = CHR_1 + ",1000,+,50M,60,2;";
        String rewritten = SaTagRewriter.rewriteSaTag(entry + entry, resolver);
        assertEquals(entry, rewritten);
    }

    @Test
    public void testAllEntriesFailReturnsNull()
    {
        LiftBackResolver resolver = newResolver();
        String sa = TX_CONTIG + ",10000,+,50M,60,0;malformed;";
        assertNull(SaTagRewriter.rewriteSaTag(sa, resolver));
    }
}
