package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.BOTH_AGREE;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.BOTH_AMBIGUOUS;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.BOTH_MULTI;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.BOTH_MULTI_TX_JUNCTION;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.REF_MULTI;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.REF_SINGLE;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.TX_MULTI;
import static com.hartwig.hmftools.tars.liftback.LiftBackCategory.TX_SINGLE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.tars.common.SpliceCommon;

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

    // The full ref-vs-tx category matrix: one assert per discriminator branch. Each assert carries a message
    // naming the branch so a failure pinpoints which category misclassified.
    @Test
    public void testCategorize()
    {
        // ---- one locus, single source ----

        // only a ref alignment present -> REF_SINGLE.
        assertEquals("single ref alignment", REF_SINGLE,
                categoryOf(set(ref(CHR1, 100, FULL_MATCH_CIGAR))));

        // only a tx (spliced) alignment present -> TX_SINGLE.
        assertEquals("single tx alignment", TX_SINGLE,
                categoryOf(set(tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));

        // ---- one locus, ref and tx both present (the discriminator) ----

        // identical cigar, no N: the two views agree, no contest -> BOTH_AGREE.
        assertEquals("ref and tx agree on the same gapless cigar", BOTH_AGREE,
                categoryOf(set(
                        ref(CHR1, 100, FULL_MATCH_CIGAR),
                        tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0))));

        // tx spliced across a real junction, ref ran off the exon edge into a softclip: tx wins.
        assertEquals("tx junction beats a ref softclip", BOTH_TX_JUNCTION_REF_SOFTCLIP,
                categoryOf(set(
                        ref(CHR1, 100, SOFTCLIP_CIGAR),
                        tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));

        // ref matched cleanly through the supposed intron (no softclip): read is unspliced, ref wins.
        assertEquals("clean ref match through the intron beats the tx junction", BOTH_TX_JUNCTION_REF_MATCH,
                categoryOf(set(
                        ref(CHR1, 100, FULL_MATCH_CIGAR),
                        tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));

        // tx softclipped at the exon boundary (no N junction), ref full match: ref wins (intron retention).
        assertEquals("ref full match beats a tx that only softclips at the boundary", BOTH_TX_SOFTCLIP_REF_MATCH,
                categoryOf(set(
                        ref(CHR1, 100, FULL_MATCH_CIGAR),
                        tx(CHR1, 100, "50M50S", true, 0))));

        // tx contiguous (no junction, no boundary clip), ref softclipped: no discriminator rule applies.
        assertEquals("no rule fires when tx is contiguous and ref is softclipped", BOTH_AMBIGUOUS,
                categoryOf(set(
                        ref(CHR1, 100, SOFTCLIP_CIGAR),
                        tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0))));

        // ---- two or more loci ----

        // multiple loci, all ref, zero tx alts -> REF_MULTI.
        assertEquals("multi-locus, ref only", REF_MULTI,
                categoryOf(set(
                        ref(CHR1, 100, FULL_MATCH_CIGAR),
                        ref(CHR2, 200, FULL_MATCH_CIGAR))));

        // multiple loci, all tx, zero ref alts -> TX_MULTI.
        assertEquals("multi-locus, tx only", TX_MULTI,
                categoryOf(set(
                        tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0),
                        tx(CHR2, 200, TX_JUNCTION_CIGAR, false, 0))));

        // distinct loci, ref and tx; tx has an annotated junction, ref alt is intronless (paralog): tx faithful.
        assertEquals("multi-locus with a tx junction against an intronless ref paralog", BOTH_MULTI_TX_JUNCTION,
                categoryOf(set(
                        tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0),
                        ref(CHR2, 200, FULL_MATCH_CIGAR))));

        // distinct loci, neither side has a junction: genuine multi-mapper.
        assertEquals("multi-locus with no junction on either side", BOTH_MULTI,
                categoryOf(set(
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
    public void testTxSoftclipRefMatchTxSelfSwapsToRef()
    {
        // self is the softclipped tx placement; a clean ref full-match sits at the same locus (intron
        // retention). Ref wins: swap primary to the ref alt, demote self.
        final LiftedAlignment self = tx(CHR1, 100, "50M50S", true, 0);
        self.IsPrimaryChoice = true;
        final LiftedAlignment refAlt = ref(CHR1, 100, FULL_MATCH_CIGAR);

        final LiftBackDiscriminator.Outcome outcome =
                LiftBackDiscriminator.apply(set(self, refAlt), BOTH_TX_SOFTCLIP_REF_MATCH, self);

        assertSame(refAlt, outcome.effectivePrimary());
        assertEquals("swapped_tx_to_ref", outcome.note());
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
