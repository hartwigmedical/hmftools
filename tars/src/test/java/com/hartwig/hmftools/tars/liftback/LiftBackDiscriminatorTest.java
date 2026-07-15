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
        for(LiftedAlignment alignment : alignments)
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

    // ---- apply(): score-based primary pick (no shape rules) ----

    // single locus: a ref (softclipped) and a tx (contiguous) candidate that only score can separate.
    private static List<LiftedAlignment> contestedSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "50M51S");
        self.IsPrimaryChoice = true;
        LiftedAlignment txAlt = tx(CHR1, 100, "100M", false, 0);
        return set(self, txAlt);
    }

    // two loci: a ref and a tx placement.
    private static List<LiftedAlignment> multiLocusSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "151M");
        self.IsPrimaryChoice = true;
        LiftedAlignment altB = tx(CHR2, 200, "151M", false, 0);
        return set(self, altB);
    }

    @Test
    public void testBwaPriorityKeepsSelf()
    {
        // MAPQ > 0: bwa ranked the placements, so its primary is kept regardless of score.
        List<LiftedAlignment> alignments = contestedSet();
        alignments.get(0).GenomicScore = 10;
        alignments.get(1).GenomicScore = 99;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, AMBIGUOUS, alignments.get(0), 0, true);
        assertSame(alignments.get(0), outcome.effectivePrimary());
        assertEquals("", outcome.note());
    }

    @Test
    public void testUnscoredKeepsSelf()
    {
        // contested but no candidate scored (a split read left for Step 3): keep bwa's primary, drop nothing.
        List<LiftedAlignment> alignments = contestedSet();
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, AMBIGUOUS, alignments.get(0), 0, false);
        assertSame(alignments.get(0), outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertFalse(alignments.get(0).Dropped);
        assertFalse(alignments.get(1).Dropped);
    }

    @Test
    public void testDecisiveScoreWinsRegardlessOfSeed()
    {
        // higher genome score wins outright; the seed does not matter.
        List<LiftedAlignment> refWins = contestedSet();
        refWins.get(0).GenomicScore = 90;
        refWins.get(1).GenomicScore = 50;
        LiftBackDiscriminator.ApplyResult refResult = LiftBackDiscriminator.apply(refWins, AMBIGUOUS, refWins.get(0), 1, false);
        assertSame(refWins.get(0), refResult.effectivePrimary());
        assertEquals("score", refResult.note());

        List<LiftedAlignment> txWins = contestedSet();
        txWins.get(0).GenomicScore = 40;
        txWins.get(1).GenomicScore = 88;
        LiftBackDiscriminator.ApplyResult txResult = LiftBackDiscriminator.apply(txWins, AMBIGUOUS, txWins.get(0), 0, false);
        assertSame(txWins.get(1), txResult.effectivePrimary());
        assertEquals("score", txResult.note());
        assertTrue(txWins.get(1).IsPrimaryChoice);
        assertFalse(txWins.get(0).IsPrimaryChoice);
        assertFalse("loser rides in XA, not dropped", txWins.get(0).Dropped);
    }

    @Test
    public void testScoreTieFallsToSeededRandom()
    {
        // equal scores -> seeded pick, reproducible. contestedSet order is [self, tx]; seed 0 -> self, seed 1 -> tx.
        List<LiftedAlignment> even = contestedSet();
        even.get(0).GenomicScore = 70;
        even.get(1).GenomicScore = 70;
        LiftBackDiscriminator.ApplyResult evenResult = LiftBackDiscriminator.apply(even, AMBIGUOUS, even.get(0), 0, false);
        assertSame(even.get(0), evenResult.effectivePrimary());
        assertEquals("random", evenResult.note());

        List<LiftedAlignment> odd = contestedSet();
        odd.get(0).GenomicScore = 70;
        odd.get(1).GenomicScore = 70;
        LiftBackDiscriminator.ApplyResult oddResult = LiftBackDiscriminator.apply(odd, AMBIGUOUS, odd.get(0), 1, false);
        assertSame(odd.get(1), oddResult.effectivePrimary());
        assertEquals("random", oddResult.note());
        assertFalse("tie loser rides in XA, not dropped", odd.get(0).Dropped);
    }

    @Test
    public void testScoreTieJunctionBeatsSoftClipAtSameLocus()
    {
        // a spliced placement and a soft-clip placement at the same locus, tied on score -> the junction wins
        // outright (not the coin), regardless of seed. set order is [soft-clip, junction] and seed 0 would pick
        // the soft-clip if this were a plain random tie, so the junction winning proves the shape rule fired.
        LiftedAlignment softClip = ref(CHR1, 100, SOFTCLIP_CIGAR);
        softClip.IsPrimaryChoice = true;
        LiftedAlignment junction = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, AMBIGUOUS, softClip, 0, false);
        assertSame(junction, outcome.effectivePrimary());
        assertEquals("junction", outcome.note());
        assertTrue(junction.IsPrimaryChoice);
        assertFalse(softClip.IsPrimaryChoice);
        assertFalse("soft-clip loser rides in XA, not dropped", softClip.Dropped);
    }

    @Test
    public void testScoreTieJunctionAtDifferentLocusStaysRandom()
    {
        // junction and soft-clip sit at different loci, so the same-locus rule does not fire and the tie is random.
        LiftedAlignment softClip = ref(CHR1, 100, SOFTCLIP_CIGAR);
        softClip.IsPrimaryChoice = true;
        LiftedAlignment junction = tx(CHR2, 200, TX_JUNCTION_CIGAR, false, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, AMBIGUOUS, softClip, 0, false);
        assertSame("seed 0 -> first candidate", softClip, outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testMultiLocusDecisiveScorePicksBestLocus()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 60;
        alignments.get(1).GenomicScore = 130;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, MULTIMAPPER, alignments.get(0), 0, false);
        assertSame(alignments.get(1), outcome.effectivePrimary());
        assertEquals("score", outcome.note());
        assertFalse("all placements ride in XA", alignments.get(0).Dropped);
    }

    @Test
    public void testMultiLocusScoreTieSeededRandom()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, MULTIMAPPER, alignments.get(0), 1, false);
        assertSame("seed 1 -> second candidate", alignments.get(1), outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testScoreTieCollapsesIdenticalPlacements()
    {
        // self and a tx alt lift to the same locus and CIGAR; a third tx alt is a distinct spliced placement. The two
        // identical placements collapse to one, so the seeded tie is over two distinct placements, not three: seed 1
        // lands on the spliced placement (it would land on the duplicate 100M without the collapse).
        LiftedAlignment self = ref(CHR1, 100, "100M");
        self.IsPrimaryChoice = true;
        LiftedAlignment txSame = tx(CHR1, 100, "100M", false, 0);
        LiftedAlignment txSpliced = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        List<LiftedAlignment> alignments = set(self, txSame, txSpliced);
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        alignments.get(2).GenomicScore = 100;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, MULTIMAPPER, self, 1, false);
        assertSame(txSpliced, outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testMultiLocusTieMateProximityPicksMateLocus()
    {
        // tied loci CHR1:100 and CHR2:200; the mate maps on CHR2 near 200, so that locus wins over the seed.
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedMateInfo mate = LiftedMateInfo.mapped(CHR2, 250, 350, "100M", false);
        // seed 0 would pick CHR1:100; mate proximity overrides to the CHR2 locus.
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, MULTIMAPPER, alignments.get(0), 0, false, mate);
        assertSame(alignments.get(1), outcome.effectivePrimary());
        assertEquals("mate", outcome.note());
        assertTrue(alignments.get(1).IsPrimaryChoice);
        assertFalse("tie loser rides in XA, not dropped", alignments.get(0).Dropped);
    }

    @Test
    public void testMultiLocusTieMateOnNeitherChromStaysRandom()
    {
        // the mate is on a third chromosome, so proximity does not discriminate and the seed decides.
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedMateInfo mate = LiftedMateInfo.mapped("chr9", 500, 600, "100M", false);
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, MULTIMAPPER, alignments.get(0), 1, false, mate);
        assertSame("seed 1 -> second candidate", alignments.get(1), outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testMultiLocusTieMateTooFarStaysRandom()
    {
        // mate is on CHR1 but past MATE_PROXIMITY_MAX_DISTANCE from the CHR1 locus, so it is not proximal.
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedMateInfo mate = LiftedMateInfo.mapped(CHR1, 5_000_000, 5_000_100, "100M", false);
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, MULTIMAPPER, alignments.get(0), 1, false, mate);
        assertSame("seed 1 -> second candidate", alignments.get(1), outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }
}
