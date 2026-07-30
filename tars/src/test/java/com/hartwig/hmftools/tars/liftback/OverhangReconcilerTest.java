package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.recordBuilder;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.selfAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

// Unit tests for OverhangReconciler: reconcileAlignmentsToGenome (per-candidate collapse, XA-alt drop, genomic re-score)
// and reconcileChosenPrimary (collapse + tx-match-only softclip extension). Gate mechanics are covered in OverhangGateTest.
public class OverhangReconcilerTest
{
    // Standard all-'A' genome, matching OverhangGateTest so the reused collapse/extend outcomes hold.
    private static TestGenome genome()
    {
        return new TestGenome().with(CHR_1, 5000, 'A');
    }

    private static LiftedAlignment refAlt(final String chrom, final int pos, final String cigar, final boolean forwardStrand)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, false, forwardStrand, 0);
    }

    // tx-contig alignment; the record builder marks index 0 as the primary, so this drives the softclip extension.
    private static LiftedAlignment txPrimaryChoice(final int pos, final String cigar)
    {
        return new LiftedAlignment(CHR_1, pos, cigar, 0, true, false, true, 1);
    }

    // -- reconcileAlignmentsToGenome ------------------------------------------------------------------------------

    @Test
    public void testSingleCandidateFastPathCollapsesSelfWithoutDropping()
    {
        // one candidate -> fast path: self's "20M100N8M4S" collapses to "20M12S" (the 8M anchor is C bases against an
        // all-'A' genome, so it scores below the gate). Self keeps the collapsed cigar and is never dropped.
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        List<LiftedAlignment> alignments = new ArrayList<>(List.of(selfAlignment(CHR_1, 1, "20M100N8M4S")));
        SAMRecord record = primaryRecord(CHR_1, 1, "20M100N8M4S");
        record.setReadBases(bases("C".repeat(32)));

        reconciler.reconcileAlignmentsToGenome(alignments, record);

        assertEquals("20M12S", alignments.get(0).LiftedCigar);
        assertEquals(1, alignments.get(0).LiftedPos);
        assertFalse(alignments.get(0).Dropped);
    }

    @Test
    public void testCollapsedXaAltIsDropped()
    {
        // the XA alt's junction fully collapses (fabricated placement) so it is flagged Dropped, but left
        // un-rewritten in the list; self (a clean ref) is never dropped.
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment alt = refAlt(CHR_1, 300, "20M100N8M4S", true);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, alt));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        record.setReadBases(bases("C".repeat(50)));

        reconciler.reconcileAlignmentsToGenome(alignments, record);

        assertTrue(alt.Dropped);
        assertEquals("20M100N8M4S", alt.LiftedCigar);
        assertFalse(self.Dropped);
    }

    @Test
    public void testGenomicScoreSetOnCandidatesWithoutSupplementary()
    {
        // no SA tag -> each surviving candidate is re-scored against the genome (50 matches vs all-'A' = 50).
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment alt = refAlt(CHR_1, 2000, "50M", true);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, alt));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        record.setReadBases(bases("A".repeat(50)));

        reconciler.reconcileAlignmentsToGenome(alignments, record);

        assertEquals(50, self.GenomicScore);
        assertEquals(50, alt.GenomicScore);
    }

    @Test
    public void testSplitReadLeavesGenomicScoreUnset()
    {
        // a same-chromosome supplementary within an intron marks a genuine split read: leave GenomicScore at the
        // sentinel so the supplementary-resolve merge stands.
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment alt = refAlt(CHR_1, 2000, "50M", true);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, alt));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        record.setReadBases(bases("A".repeat(50)));
        record.setAttribute("SA", CHR_1 + ",2000,+,50M,0,0;");

        reconciler.reconcileAlignmentsToGenome(alignments, record);

        assertEquals(Integer.MIN_VALUE, self.GenomicScore);
        assertEquals(Integer.MIN_VALUE, alt.GenomicScore);
    }

    @Test
    public void testCrossChromSupplementaryIsScored()
    {
        // a cross-chromosome supplementary is never merged, so the contiguous primary IS scored.
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment alt = refAlt(CHR_1, 2000, "50M", true);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, alt));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        record.setReadBases(bases("A".repeat(50)));
        record.setAttribute("SA", "chr2,500,+,20M30S,0,0;");   // supplementary on a different chromosome

        reconciler.reconcileAlignmentsToGenome(alignments, record);

        assertEquals(50, self.GenomicScore);
        assertEquals(50, alt.GenomicScore);
    }

    @Test
    public void testFarSameChromSupplementaryIsScored()
    {
        // a same-chromosome supplementary beyond the max intron span is not a splice partner, so the primary is scored.
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment alt = refAlt(CHR_1, 2000, "50M", true);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, alt));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        record.setReadBases(bases("A".repeat(50)));
        record.setAttribute("SA", CHR_1 + ",50000000,+,20M30S,0,0;");   // same chrom but 50 Mb away (> max intron)

        reconciler.reconcileAlignmentsToGenome(alignments, record);

        assertEquals(50, self.GenomicScore);
        assertEquals(50, alt.GenomicScore);
    }

    @Test
    public void testOppositeStrandAltScoredOnReverseComplement()
    {
        // the alt is opposite-strand to the forward record, so the gate scores it on the reverse-complemented read.
        // The all-'A' read reverse-complements to all-'T', matching the all-'T' patch at chr1:2000; the forward bases
        // would mismatch it. GenomicScore therefore equals the reverse-complement score, and beats the forward one.
        TestGenome genome = genome().set(CHR_1, 2000, "T".repeat(50));
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment alt = refAlt(CHR_1, 2000, "50M", false);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, alt));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        byte[] forward = bases("A".repeat(50));
        record.setReadBases(forward);

        reconciler.reconcileAlignmentsToGenome(alignments, record);

        byte[] reverseComplement = Arrays.copyOf(forward, forward.length);
        SequenceUtil.reverseComplement(reverseComplement);
        OverhangGate scorer = new OverhangGate(ref);
        int reverseScore = scorer.genomicScore(CHR_1, 2000, "50M", reverseComplement);
        int forwardScore = scorer.genomicScore(CHR_1, 2000, "50M", forward);

        assertEquals(reverseScore, alt.GenomicScore);
        assertTrue(forwardScore < reverseScore);
    }

    // -- reconcileChosenPrimary -----------------------------------------------------------------------------------

    @Test
    public void testCollapsesFabricatedJunctionOnChosenPrimary()
    {
        // the chosen (tx) primary's fabricated "20M100N8M4S" merge collapses to "20M12S" - the collapse catches a junction
        // fabricated after resolve, not just the softclip extension.
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedRecord primaryRes = recordBuilder()
                .alignments(new ArrayList<>(List.of(txPrimaryChoice(1, "20M100N8M4S"))))
                .build();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "20M100N8M4S");
        primary.setReadBases(bases("C".repeat(32)));

        LiftedRecord[] resolved = { primaryRes };
        reconciler.reconcileChosenPrimary(resolved, primary, 0, null);

        assertEquals("20M12S", resolved[0].finalCigar());
        assertEquals(1, resolved[0].finalPos());
    }

    // Ultima single-end: the reconcile path reads strand and bases but never a pairing flag, so an unpaired record
    // extends its tail exactly as the paired case does.
    @Test
    public void testSingleEndPrimaryExtendsTerminalSoftClip()
    {
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedRecord primaryRes = recordBuilder()
                .alignments(new ArrayList<>(List.of(txPrimaryChoice(1, "50M10S"))))
                .build();
        SAMRecord primary = TarsTestFixtures.unpairedPrimaryRecord(TX_CONTIG, 1, "50M10S");
        primary.setReadBases(bases("A".repeat(60)));

        LiftedRecord[] resolved = { primaryRes };
        reconciler.reconcileChosenPrimary(resolved, primary, 0, null);

        assertEquals("60M", resolved[0].finalCigar());
    }

    @Test
    public void testTxMatchPrimaryExtendsTerminalSoftClip()
    {
        // tx-match primary with a "50M10S" tail whose clip continues contiguously in the genome (all-'A') extends the
        // softclip to "60M".
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedRecord primaryRes = recordBuilder()
                .alignments(new ArrayList<>(List.of(txPrimaryChoice(1, "50M10S"))))
                .build();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M10S");
        primary.setReadBases(bases("A".repeat(60)));

        LiftedRecord[] resolved = { primaryRes };
        reconciler.reconcileChosenPrimary(resolved, primary, 0, null);

        assertEquals("60M", resolved[0].finalCigar());
        assertEquals(1, resolved[0].finalPos());
    }

    @Test
    public void testRefPrimaryTerminalSoftClipNotExtended()
    {
        // identical "50M10S" tail and genome as the tx case, but the chosen primary is a genomic/ref placement: the
        // the extension is gated off, so bwa's clip is left untouched and the result is not revised.
        TestGenome genome = genome();
        RefGenomeInterface ref = genome.asRefGenome();
        OverhangReconciler reconciler = new OverhangReconciler(new OverhangGate(ref));

        LiftedAlignment refPrimary = selfAlignment(CHR_1, 1, "50M10S");
        LiftedRecord primaryRes = recordBuilder()
                .alignments(new ArrayList<>(List.of(refPrimary)))
                .build();
        SAMRecord primary = primaryRecord(CHR_1, 1, "50M10S");
        primary.setReadBases(bases("A".repeat(60)));

        LiftedRecord[] resolved = { primaryRes };
        reconciler.reconcileChosenPrimary(resolved, primary, 0, null);

        assertEquals("50M10S", resolved[0].finalCigar());
        assertSame(primaryRes, resolved[0]);
    }

    @Test
    public void testNullGateLeavesChosenPrimaryUntouched()
    {
        // no gate loaded: reconcileChosenPrimary is a no-op even for a tx primary whose "50M10S" tail would extend.
        OverhangReconciler reconciler = new OverhangReconciler(null);

        LiftedRecord primaryRes = recordBuilder()
                .alignments(new ArrayList<>(List.of(txPrimaryChoice(1, "50M10S"))))
                .build();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M10S");
        primary.setReadBases(bases("A".repeat(60)));

        LiftedRecord[] resolved = { primaryRes };
        reconciler.reconcileChosenPrimary(resolved, primary, 0, null);

        assertSame(primaryRes, resolved[0]);
        assertEquals("50M10S", resolved[0].finalCigar());
    }

    @Test
    public void testDisabledGateLeavesChosenPrimaryUntouched()
    {
        // gate present but disabled (no ref handle): the enabled() guard skips the collapse/extend block, so the same
        // extendable tx primary is left as-is.
        OverhangGate disabledGate = new OverhangGate(null);
        OverhangReconciler reconciler = new OverhangReconciler(disabledGate);

        LiftedRecord primaryRes = recordBuilder()
                .alignments(new ArrayList<>(List.of(txPrimaryChoice(1, "50M10S"))))
                .build();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M10S");
        primary.setReadBases(bases("A".repeat(60)));

        LiftedRecord[] resolved = { primaryRes };
        reconciler.reconcileChosenPrimary(resolved, primary, 0, null);

        assertSame(primaryRes, resolved[0]);
        assertEquals("50M10S", resolved[0].finalCigar());
    }
}
