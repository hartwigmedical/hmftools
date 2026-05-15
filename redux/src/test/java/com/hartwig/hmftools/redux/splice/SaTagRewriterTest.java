package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

// covers SaTagRewriter: ref pass-through, Tx-contig lift, malformed entry skip,
// dedup of post-lift duplicates, null/empty input, and the all-fail -> null contract.
public class SaTagRewriterTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    // same three-exon contig as LiftBackResolverTest: exons chr1:100-199, 300-399, 500-549; contig length 250.
    private static LiftBackResolver newResolver()
    {
        ContigEntry entry = new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1,
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
        // alt-contig pos 1 corresponds to exon 1 start (chr1:100). 50M stays 50M (no junction crossed).
        String sa = TX_CONTIG + ",1,+,50M,60,2;";
        String rewritten = SaTagRewriter.rewriteSaTag(sa, resolver);
        assertEquals(CHR_1 + ",100,+,50M,60,2;", rewritten);
    }

    @Test
    public void testMalformedEntriesSkipped()
    {
        LiftBackResolver resolver = newResolver();
        // too few fields, non-numeric position, then a valid ref entry
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
        // position past exon 3 end on the tx contig -> outside any transcript span -> findSegment returns null
        String sa = TX_CONTIG + ",10000,+,50M,60,0;malformed;";
        assertNull(SaTagRewriter.rewriteSaTag(sa, resolver));
    }
}
