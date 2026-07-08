package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Set;

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
        LiftBackResolver resolver = newResolver();

        // null or empty input -> null
        assertNull(SaTagRewriter.rewriteSaTag(null, resolver));
        assertNull(SaTagRewriter.rewriteSaTag("", resolver));

        // ref-contig entry passes through unchanged
        String refEntry = CHR_1 + ",1000,+,50M,60,2;";
        assertEquals(refEntry, SaTagRewriter.rewriteSaTag(refEntry, resolver));

        // tx-contig entry lifted to genomic coordinates
        assertEquals(CHR_1 + ",100,+,50M,60,2;", SaTagRewriter.rewriteSaTag(TX_CONTIG + ",1,+,50M,60,2;", resolver));

        // a malformed entry drops the whole tag (parsing is delegated to hmf-common's SupplementaryReadData,
        // which rejects the tag rather than skipping bad entries; bwa never emits malformed SA in practice)
        assertNull(SaTagRewriter.rewriteSaTag("junk;" + refEntry, resolver));

        // duplicate lifted entries deduped
        assertEquals(refEntry, SaTagRewriter.rewriteSaTag(refEntry + refEntry, resolver));

        // all entries fail -> null
        assertNull(SaTagRewriter.rewriteSaTag(TX_CONTIG + ",10000,+,50M,60,0;malformed;", resolver));
    }

    @Test
    public void testExcludeKeyRemovesDroppedSuppAcrossClipType()
    {
        LiftBackResolver resolver = newResolver();

        // The primary's SA lists a dropped supp soft-clipped (19M247S); the exclude key (built by suppSaKey from the
        // hard-clipped supp record, normalised H->S) is soft. The entry must be removed, not left dangling.
        Set<String> exclude = Set.of(CHR_1 + ":1000:-:19M247S");
        assertNull(SaTagRewriter.rewriteSaTag(CHR_1 + ",1000,-,19M247S,0,0;", resolver, exclude));

        // and when the SA entry itself is hard-clipped, rewriteSaTag normalises H->S before matching the same key.
        assertNull(SaTagRewriter.rewriteSaTag(CHR_1 + ",1000,-,19M247H,0,0;", resolver, exclude));

        // an entry at a different locus is untouched by the exclude set
        assertEquals(CHR_1 + ",2000,-,19M247S,0,0;",
                SaTagRewriter.rewriteSaTag(CHR_1 + ",2000,-,19M247S,0,0;", resolver, exclude));
    }
}
