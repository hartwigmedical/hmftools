package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.DecidingFeature.AMBIGUOUS;
import static com.hartwig.hmftools.tars.liftback.DecidingFeature.CONCORDANT;
import static com.hartwig.hmftools.tars.liftback.DecidingFeature.JUNCTION;
import static com.hartwig.hmftools.tars.liftback.DecidingFeature.JUNCTION_OVER_CONTIGUOUS;
import static com.hartwig.hmftools.tars.liftback.DecidingFeature.MULTIMAPPER;
import static com.hartwig.hmftools.tars.liftback.DecidingFeature.REF_READS_THROUGH;
import static com.hartwig.hmftools.tars.liftback.DecidingFeature.SOLE_REF;
import static com.hartwig.hmftools.tars.liftback.DecidingFeature.SOLE_TX;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

// Direct unit coverage of the ref-vs-tx discriminator: categorize() resolves a record's lifted alignment set
// to an Outcome + DecidingFeature, and apply() executes the swap/drop that feature implies. These are the
// central manipulations that decide whether a tx-contig (spliced) placement wins over bwa's ref placement.
// Each test feeds a hand-built alignment set so the decision and its consequence are asserted in isolation,
// rather than through the full LiftBackResolver.resolve path.
public class LiftBackDiscriminatorTest
{
    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";

    // anchors comfortably above TarsConstants.MIN_JUNCTION_ANCHOR so the N counts as a real junction.
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
        List<LiftedAlignment> list = new ArrayList<>();
        for(final LiftedAlignment alignment : alignments)
        {
            list.add(alignment);
        }
        return list;
    }

    private static DecidingFeature featureOf(final List<LiftedAlignment> alignments)
    {
        return LiftBackDiscriminator.categorize(alignments).Feature;
    }

    // The full ref-vs-tx decision matrix: one assert per discriminator branch. Each assert carries a message
    // naming the branch so a failure pinpoints which decision misclassified.
    @Test
    public void testCategorize()
    {
        // ---- one locus, single source ----

        // only a ref alignment present -> REF, SOLE_REF.
        assertEquals("single ref alignment", SOLE_REF,
                featureOf(set(ref(CHR1, 100, FULL_MATCH_CIGAR))));

        // only a tx (spliced) alignment present -> TX, SOLE_TX.
        assertEquals("single tx alignment", SOLE_TX,
                featureOf(set(tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));

        // ---- one locus, ref and tx both present (the discriminator) ----

        // identical cigar, no N: the two views agree, no contest -> CONCORDANT.
        assertEquals("ref and tx agree on the same gapless cigar", CONCORDANT,
                featureOf(set(
                        ref(CHR1, 100, FULL_MATCH_CIGAR),
                        tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0))));

        // tx spliced across a real junction, ref ran off the exon edge into a softclip: tx wins (JUNCTION).
        assertEquals("tx junction beats a ref softclip", JUNCTION,
                featureOf(set(
                        ref(CHR1, 100, SOFTCLIP_CIGAR),
                        tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0))));

        // ref matched cleanly through the supposed intron (no softclip): read is unspliced, ref wins.
        // REF_READS_THROUGH merges the old TX_JUNCTION_REF_MATCH / TX_SOFTCLIP_REF_MATCH; the underlying
        // flags keep the two scenarios distinguishable, so assert them too.
        LiftBackDiscriminator.Features refMatch = LiftBackDiscriminator.categorize(set(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0)));
        assertEquals("clean ref match through the intron beats the tx junction", REF_READS_THROUGH, refMatch.Feature);
        assertTrue(refMatch.TxHasNCigar);
        assertTrue(refMatch.RefFullMatch);
        assertFalse(refMatch.RefSoftClipped);

        // tx softclipped at the exon boundary (no N junction), ref full match: ref wins (intron retention).
        LiftBackDiscriminator.Features txSoftclip = LiftBackDiscriminator.categorize(set(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, "50M50S", true, 0)));
        assertEquals("ref full match beats a tx that only softclips at the boundary", REF_READS_THROUGH, txSoftclip.Feature);
        assertTrue(txSoftclip.TxSoftClipAtBoundary);
        assertFalse(txSoftclip.TxHasNCigar);

        // tx contiguous (no junction, no boundary clip), ref softclipped: no rule applies -> AMBIGUOUS.
        assertEquals("no rule fires when tx is contiguous and ref is softclipped", AMBIGUOUS,
                featureOf(set(
                        ref(CHR1, 100, SOFTCLIP_CIGAR),
                        tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0))));

        // ---- two or more loci ----

        // multiple loci, all ref, zero tx alts -> SOLE_REF.
        assertEquals("multi-locus, ref only", SOLE_REF,
                featureOf(set(
                        ref(CHR1, 100, FULL_MATCH_CIGAR),
                        ref(CHR2, 200, FULL_MATCH_CIGAR))));

        // multiple loci, all tx, zero ref alts -> SOLE_TX.
        assertEquals("multi-locus, tx only", SOLE_TX,
                featureOf(set(
                        tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0),
                        tx(CHR2, 200, TX_JUNCTION_CIGAR, false, 0))));

        // distinct loci, ref and tx; tx has an annotated junction, ref alt is intronless (intronless ref alt): tx faithful.
        assertEquals("multi-locus with a tx junction against an intronless ref alt", JUNCTION_OVER_CONTIGUOUS,
                featureOf(set(
                        tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0),
                        ref(CHR2, 200, FULL_MATCH_CIGAR))));

        // distinct loci, neither side has a junction: genuine multi-mapper.
        assertEquals("multi-locus with no junction on either side", MULTIMAPPER,
                featureOf(set(
                        tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0),
                        ref(CHR2, 200, FULL_MATCH_CIGAR))));
    }

    @Test
    public void testCategorizeSkipsGateDroppedAlt()
    {
        // an XA alt the overhang gate collapsed to a contiguous alignment is marked Dropped before the
        // discriminator runs; categorize must ignore it. A ref self plus a dropped tx alt at another locus
        // resolves to SOLE_REF (one counted locus), not a multi-locus tx-vs-ref contest.
        LiftedAlignment self = ref(CHR1, 100, FULL_MATCH_CIGAR);
        LiftedAlignment droppedTx = tx(CHR2, 200, TX_JUNCTION_CIGAR, false, 0);
        droppedTx.Dropped = true;

        assertEquals(SOLE_REF, featureOf(set(self, droppedTx)));
    }

    // ---- apply(): tx-favouring swap (promoteTxOverRef) ----

    @Test
    public void testTxJunctionRefSoftclipSwapsRefSelfToTx()
    {
        // bwa picked the ref placement; tx wins, so the primary swaps to the tx alt. self stays as an alt
        // (preserves its locus) but loses IsPrimaryChoice; any other ref alt is dropped.
        LiftedAlignment self = ref(CHR1, 100, SOFTCLIP_CIGAR);
        self.IsPrimaryChoice = true;
        LiftedAlignment txWinner = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        LiftedAlignment otherRef = ref(CHR2, 200, FULL_MATCH_CIGAR);

        LiftBackDiscriminator.ApplyResult outcome =
                LiftBackDiscriminator.apply(set(self, txWinner, otherRef), JUNCTION, self);

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
        LiftedAlignment self = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        self.IsPrimaryChoice = true;
        LiftedAlignment refAlt = ref(CHR1, 100, SOFTCLIP_CIGAR);

        LiftBackDiscriminator.ApplyResult outcome =
                LiftBackDiscriminator.apply(set(self, refAlt), JUNCTION, self);

        assertSame(self, outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertTrue(refAlt.Dropped);
    }

    @Test
    public void testTxSwapPrefersFewestMismatchJunctionAlt()
    {
        // two tx junction alts at the locus: the winner is the N-junction alt with the fewest mismatches.
        LiftedAlignment self = ref(CHR1, 100, SOFTCLIP_CIGAR);
        self.IsPrimaryChoice = true;
        LiftedAlignment txMany = tx(CHR1, 100, "50M100N50M", false, 5);
        LiftedAlignment txFew = tx(CHR1, 100, "60M100N40M", false, 2);

        LiftBackDiscriminator.ApplyResult outcome =
                LiftBackDiscriminator.apply(set(self, txMany, txFew), JUNCTION, self);

        assertSame(txFew, outcome.effectivePrimary());
    }

    // ---- apply(): ref-favouring swap (promoteRefOverTx) ----

    @Test
    public void testTxJunctionRefMatchRefSelfDropsTx()
    {
        // ref matched through the intron; self is ref, so keep it and drop the tx alt.
        LiftedAlignment self = ref(CHR1, 100, FULL_MATCH_CIGAR);
        self.IsPrimaryChoice = true;
        LiftedAlignment txAlt = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);

        LiftBackDiscriminator.ApplyResult outcome =
                LiftBackDiscriminator.apply(set(self, txAlt), REF_READS_THROUGH, self);

        assertSame(self, outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertTrue(txAlt.Dropped);
    }

    @Test
    public void testTxJunctionRefMatchTxSelfSwapsToRef()
    {
        // self is the tx placement but ref wins: swap the primary to the ref alt.
        LiftedAlignment self = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        self.IsPrimaryChoice = true;
        LiftedAlignment refWinner = ref(CHR1, 100, FULL_MATCH_CIGAR);

        LiftBackDiscriminator.ApplyResult outcome =
                LiftBackDiscriminator.apply(set(self, refWinner), REF_READS_THROUGH, self);

        assertSame(refWinner, outcome.effectivePrimary());
        assertEquals("swapped_tx_to_ref", outcome.note());
        assertTrue(refWinner.IsPrimaryChoice);
        assertFalse(self.IsPrimaryChoice);
    }

    // ---- apply(): tx-softclip-ref-match (also REF_READS_THROUGH) ----

    @Test
    public void testTxSoftclipRefMatchRefSelfDropsTx()
    {
        LiftedAlignment self = ref(CHR1, 100, FULL_MATCH_CIGAR);
        self.IsPrimaryChoice = true;
        LiftedAlignment txAlt = tx(CHR1, 100, "50M50S", true, 0);

        LiftBackDiscriminator.ApplyResult outcome =
                LiftBackDiscriminator.apply(set(self, txAlt), REF_READS_THROUGH, self);

        assertSame(self, outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertTrue(txAlt.Dropped);
    }

    @Test
    public void testTxSoftclipRefMatchTxSelfSwapsToRef()
    {
        // self is the softclipped tx placement; a clean ref full-match sits at the same locus (intron
        // retention). Ref wins: swap primary to the ref alt, demote self.
        LiftedAlignment self = tx(CHR1, 100, "50M50S", true, 0);
        self.IsPrimaryChoice = true;
        LiftedAlignment refAlt = ref(CHR1, 100, FULL_MATCH_CIGAR);

        LiftBackDiscriminator.ApplyResult outcome =
                LiftBackDiscriminator.apply(set(self, refAlt), REF_READS_THROUGH, self);

        assertSame(refAlt, outcome.effectivePrimary());
        assertEquals("swapped_tx_to_ref", outcome.note());
        assertFalse(refAlt.Dropped);
    }

    // ---- apply(): features with no swap ----

    @Test
    public void testNonSwappingFeaturesLeaveSelfAndDropNothing()
    {
        for(final DecidingFeature feature : new DecidingFeature[] {
                SOLE_REF, SOLE_TX, CONCORDANT, AMBIGUOUS, MULTIMAPPER })
        {
            LiftedAlignment self = ref(CHR1, 100, FULL_MATCH_CIGAR);
            self.IsPrimaryChoice = true;
            LiftedAlignment txAlt = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);

            LiftBackDiscriminator.ApplyResult outcome =
                    LiftBackDiscriminator.apply(set(self, txAlt), feature, self);

            assertSame("feature " + feature, self, outcome.effectivePrimary());
            assertEquals("feature " + feature, "", outcome.note());
            assertFalse("feature " + feature, self.Dropped);
            assertFalse("feature " + feature, txAlt.Dropped);
        }
    }
}
