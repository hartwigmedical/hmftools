package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.splice.LiftBackCategory.BOTH_AGREE;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.BOTH_AMBIGUOUS;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.BOTH_MULTI;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.BOTH_MULTI_TX_JUNCTION;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.REF_MULTI;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.REF_SINGLE;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.TX_MULTI;
import static com.hartwig.hmftools.redux.splice.LiftBackCategory.TX_SINGLE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

// Direct unit coverage of the ref-vs-tx discriminator: categorize() assigns one LiftBackCategory to a
// record's lifted alignment set, and apply() executes the swap/drop that category implies. These are the
// central manipulations that decide whether a tx-contig (spliced) placement wins over bwa's ref placement.
// Each test feeds a hand-built alignment set so the category and its consequence are asserted in isolation,
// rather than through the full LiftBackResolver.resolve path.
public class LiftBackDiscriminatorTest
{
    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";

    // anchors comfortably above SpliceCommon.MIN_JUNCTION_ANCHOR so the N counts as a real junction.
    private static final String TX_JUNCTION_CIGAR = "50M100N50M";
    private static final String FULL_MATCH_CIGAR = "100M";
    private static final String SOFTCLIP_CIGAR = "50M51S";

    private static LiftedAlignment tx(
            final String chrom, final int pos, final String cigar, final boolean softClipAtBoundary, final int numMismatches)
    {
        return new LiftedAlignment(
                LiftedAlignment.AlignmentSource.XA_INPUT, "ctg_tx", pos, cigar,
                chrom, pos, cigar, 0, numMismatches,
                "ENST001", "ENSG001", "GENE", softClipAtBoundary, true);
    }

    private static LiftedAlignment ref(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, chrom, pos, cigar,
                chrom, pos, cigar, 0, 0,
                null, null, null, false, true);
    }

    private static List<LiftedAlignment> set(final LiftedAlignment... alignments)
    {
        final List<LiftedAlignment> list = new ArrayList<>();
        for(final LiftedAlignment alignment : alignments)
            list.add(alignment);
        return list;
    }

    private static LiftBackCategory categoryOf(final List<LiftedAlignment> alignments)
    {
        return LiftBackDiscriminator.categorize(alignments).Category;
    }

    // ---- categorize(): one locus, single source ----

    @Test
    public void testRefSingle()
    {
        assertEquals(REF_SINGLE, categoryOf(set(ref(CHR1, 100, FULL_MATCH_CIGAR))));
    }

    @Test
    public void testTxSingle()
    {
        assertEquals(TX_SINGLE, categoryOf(set(tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));
    }

    // ---- categorize(): one locus, ref and tx both present ----

    @Test
    public void testBothAgreeSameCigarNoJunction()
    {
        // identical cigar, no N — the two views agree, no contest.
        assertEquals(BOTH_AGREE, categoryOf(set(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0))));
    }

    @Test
    public void testBothTxJunctionRefSoftclip()
    {
        // tx spliced across a real junction, ref ran off the exon edge into a softclip -> tx wins.
        assertEquals(BOTH_TX_JUNCTION_REF_SOFTCLIP, categoryOf(set(
                ref(CHR1, 100, SOFTCLIP_CIGAR),
                tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));
    }

    @Test
    public void testBothTxJunctionRefMatch()
    {
        // ref matched cleanly through the supposed intron (no softclip) -> read is unspliced, ref wins.
        assertEquals(BOTH_TX_JUNCTION_REF_MATCH, categoryOf(set(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));
    }

    @Test
    public void testBothTxSoftclipRefMatch()
    {
        // tx softclipped at the exon boundary (no N junction), ref full match -> ref wins (intron retention).
        assertEquals(BOTH_TX_SOFTCLIP_REF_MATCH, categoryOf(set(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, "50M50S", true, 0))));
    }

    @Test
    public void testBothAmbiguousWhenNoRuleFires()
    {
        // tx contiguous (no junction, no boundary clip), ref softclipped -> no discriminator rule applies.
        assertEquals(BOTH_AMBIGUOUS, categoryOf(set(
                ref(CHR1, 100, SOFTCLIP_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0))));
    }

    // ---- categorize(): two or more loci ----

    @Test
    public void testRefMulti()
    {
        assertEquals(REF_MULTI, categoryOf(set(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                ref(CHR2, 200, FULL_MATCH_CIGAR))));
    }

    @Test
    public void testTxMulti()
    {
        assertEquals(TX_MULTI, categoryOf(set(
                tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0),
                tx(CHR2, 200, TX_JUNCTION_CIGAR, false, 0))));
    }

    @Test
    public void testBothMultiTxJunction()
    {
        // distinct loci; tx has an annotated junction, the ref alt is intronless (paralog/pseudogene) -> tx faithful.
        assertEquals(BOTH_MULTI_TX_JUNCTION, categoryOf(set(
                tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0),
                ref(CHR2, 200, FULL_MATCH_CIGAR))));
    }

    @Test
    public void testBothMultiWhenNoTxJunction()
    {
        // distinct loci, neither side has a junction -> genuine multi-mapper.
        assertEquals(BOTH_MULTI, categoryOf(set(
                tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0),
                ref(CHR2, 200, FULL_MATCH_CIGAR))));
    }

    // ---- apply(): tx-favouring swap (promoteTxOverRef) ----

    @Test
    public void testTxJunctionRefSoftclipSwapsRefSelfToTx()
    {
        // bwa picked the ref placement; tx wins, so the primary swaps to the tx alt. self stays as an alt
        // (preserves its locus) but loses IsPrimaryChoice; any other ref alt is dropped.
        final LiftedAlignment self = ref(CHR1, 100, SOFTCLIP_CIGAR);
        self.IsPrimaryChoice = true;
        final LiftedAlignment txWinner = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        final LiftedAlignment otherRef = ref(CHR2, 200, FULL_MATCH_CIGAR);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, txWinner, otherRef), BOTH_TX_JUNCTION_REF_SOFTCLIP, self);

        assertSame(txWinner, outcome.effectivePrimary());
        assertEquals("swapped_ref_to_tx", outcome.note());
        assertTrue(txWinner.IsPrimaryChoice);
        assertFalse(self.IsPrimaryChoice);
        assertFalse(self.Dropped);
        assertTrue(otherRef.Dropped);
    }

    @Test
    public void testTxJunctionRefSoftclipTxSelfDropsRef()
    {
        // self is already the tx placement; keep it and drop the ref alt.
        final LiftedAlignment self = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        self.IsPrimaryChoice = true;
        final LiftedAlignment refAlt = ref(CHR1, 100, SOFTCLIP_CIGAR);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, refAlt), BOTH_TX_JUNCTION_REF_SOFTCLIP, self);

        assertSame(self, outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertTrue(refAlt.Dropped);
    }

    @Test
    public void testTxSwapPrefersFewestMismatchJunctionAlt()
    {
        // two tx junction alts at the locus: the winner is the N-junction alt with the fewest mismatches.
        final LiftedAlignment self = ref(CHR1, 100, SOFTCLIP_CIGAR);
        self.IsPrimaryChoice = true;
        final LiftedAlignment txMany = tx(CHR1, 100, "50M100N50M", false, 5);
        final LiftedAlignment txFew = tx(CHR1, 100, "60M100N40M", false, 2);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, txMany, txFew), BOTH_TX_JUNCTION_REF_SOFTCLIP, self);

        assertSame(txFew, outcome.effectivePrimary());
    }

    // ---- apply(): ref-favouring swap (promoteRefOverTx) ----

    @Test
    public void testTxJunctionRefMatchRefSelfDropsTx()
    {
        // ref matched through the intron; self is ref, so keep it and drop the tx alt.
        final LiftedAlignment self = ref(CHR1, 100, FULL_MATCH_CIGAR);
        self.IsPrimaryChoice = true;
        final LiftedAlignment txAlt = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, txAlt), BOTH_TX_JUNCTION_REF_MATCH, self);

        assertSame(self, outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertTrue(txAlt.Dropped);
    }

    @Test
    public void testTxJunctionRefMatchTxSelfSwapsToRef()
    {
        // self is the tx placement but ref wins: swap the primary to the ref alt.
        final LiftedAlignment self = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        self.IsPrimaryChoice = true;
        final LiftedAlignment refWinner = ref(CHR1, 100, FULL_MATCH_CIGAR);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, refWinner), BOTH_TX_JUNCTION_REF_MATCH, self);

        assertSame(refWinner, outcome.effectivePrimary());
        assertEquals("swapped_tx_to_ref", outcome.note());
        assertTrue(refWinner.IsPrimaryChoice);
        assertFalse(self.IsPrimaryChoice);
    }

    // ---- apply(): tx-softclip-ref-match ----

    @Test
    public void testTxSoftclipRefMatchRefSelfDropsTx()
    {
        final LiftedAlignment self = ref(CHR1, 100, FULL_MATCH_CIGAR);
        self.IsPrimaryChoice = true;
        final LiftedAlignment txAlt = tx(CHR1, 100, "50M50S", true, 0);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, txAlt), BOTH_TX_SOFTCLIP_REF_MATCH, self);

        assertSame(self, outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertTrue(txAlt.Dropped);
    }

    @Test
    public void testTxSoftclipRefMatchTxSelfNoSwap()
    {
        // self is tx and ref-match wins, but there is nothing to promote to (self is the only mapping kept):
        // no swap, flagged for diagnostics.
        final LiftedAlignment self = tx(CHR1, 100, "50M50S", true, 0);
        self.IsPrimaryChoice = true;
        final LiftedAlignment refAlt = ref(CHR1, 100, FULL_MATCH_CIGAR);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, refAlt), BOTH_TX_SOFTCLIP_REF_MATCH, self);

        assertSame(self, outcome.effectivePrimary());
        assertEquals("self_was_tx_no_swap", outcome.note());
        assertFalse(refAlt.Dropped);
    }

    // ---- apply(): categories with no swap ----

    @Test
    public void testNonSwappingCategoriesLeaveSelfAndDropNothing()
    {
        for(final LiftBackCategory category : new LiftBackCategory[] {
                REF_SINGLE, TX_SINGLE, BOTH_AGREE, BOTH_AMBIGUOUS, REF_MULTI, TX_MULTI, BOTH_MULTI })
        {
            final LiftedAlignment self = ref(CHR1, 100, FULL_MATCH_CIGAR);
            self.IsPrimaryChoice = true;
            final LiftedAlignment txAlt = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);

            final LiftBackDiscriminator.Outcome outcome =
                    LiftBackDiscriminator.apply(set(self, txAlt), category, self);

            assertSame("category " + category, self, outcome.effectivePrimary());
            assertEquals("category " + category, "", outcome.note());
            assertFalse("category " + category, self.Dropped);
            assertFalse("category " + category, txAlt.Dropped);
        }
    }
}
