package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.exonRegionIndex;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.recordBuilder;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// Covers lifting, the primary pick and the MAPQ policy. Cigar translation itself is covered in ContigTranslatorTest.
public class PlacementSelectorTest
{
    private static List<ContigEntry> contigMap()
    {
        return List.of(threeExonContig());
    }

    // No pairing flags, and the read name seeds the random tie-break.
    private static SAMRecord newRecord(final String contig, final int pos, final String cigar)
    {
        return TarsTestFixtures.mappedRecord("read", contig, pos, cigar);
    }

    private static SAMRecord newUnmappedRecord()
    {
        return TarsTestFixtures.unmappedRecord("read");
    }

    // Ultima single-end: with no pairing flags, an unguarded htsjdk pair getter on the resolve path throws.
    // The lifted placement must match the paired case: pairing plays no part in the decision.
    @Test
    public void testSingleEndPrimaryLifted()
    {
        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(TarsTestFixtures.unpairedPrimaryRecord(TX_CONTIG, 51, "100M"));

        assertTrue(result.hasPlacement());
        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
    }

    @Test
    public void testSingleEndSupplementaryLifted()
    {
        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(
                TarsTestFixtures.unpairedSupplementaryRecord(TX_CONTIG, 51, "100M", TX_CONTIG + ",1,+,60M40S,60,0;"));

        assertTrue(result.hasPlacement());
        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
    }

    @Test
    public void testUnmapped()
    {
        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(newUnmappedRecord());

        assertFalse(result.hasPlacement());
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testRefOnlyExonic()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "150M");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(1000, result.finalPos());
        assertEquals("150M", result.finalCigar());
        assertFalse(result.hasNCigar());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxPrimaryUniqueWithinExon()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // contig pos 1, 50M -> exon 1 (chr1:100-149)

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(100, result.finalPos());
        assertEquals("50M", result.finalCigar());
        assertFalse(result.hasNCigar());
    }

    @Test
    public void testTxPrimaryUniqueJunctionCrosser()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M"); // crosses exon 1 -> exon 2 (intron 200-299)

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
    }

    @Test
    public void testRefTxAgree()
    {
        SAMRecord record = newRecord(CHR_1, 100, "50M"); // Tx alt at contig pos 1 lifts to the same locus and CIGAR
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testIntronRetRefBetter()
    {
        // ref chr1:170 50M has its last 20bp in the intron; Tx pos 71 30M20S lifts soft-clipped at the exon boundary.
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,30M20S,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxSoftClipNotAtBoundaryFallsToAmbiguous()
    {
        // Tx 25M ends mid-exon (chr1:194), so the trailing clip is not at an exon boundary -> AMBIGUOUS.
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,25M25S,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testMultiLocusTwoLoci()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // Tx primary + ref alt on a different chrom -> two loci
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.numLoci());
    }

    @Test
    public void testDistinctLocusAltBlocksMapQualityRescue()
    {
        // Two distinct genomic loci keep the read a multimapper even though the XA alt scores worse: TARS does not
        // override bwa's MAPQ 0 with a weaker reconstructed score.
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setMappingQuality(0);
        record.setAttribute("XA", "chr5,+5000,50M,3;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.numLoci());
        assertEquals(0, result.updatedMapQuality());
    }

    @Test
    public void testNestedSpanAltCountsAsSingleLocusAndBumpsMapQuality()
    {
        // A 5'-softclipped isoform alt starts at a downstream exon of the same placement: different start, but a span
        // nested inside the junction-crossing primary's, so it is one locus and the MAPQ-0 read bumps to 60.
        // Keying loci on exact start counted this as two loci and held MAPQ at 0.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M"); // -> chr1:150 50M100N50M, genomic span 150-349
        record.setMappingQuality(0);
        record.setAttribute("XA", TX_CONTIG + ",+101,50S50M,0;"); // -> chr1:300 50S50M, span 300-349 (nested)

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
        assertEquals(60, result.updatedMapQuality());
    }

    @Test
    public void testChainedOverlapAltNotMergedThroughPrimary()
    {
        // Primary 1000-1099; alt B 1080-1179 overlaps it; alt C 1160-1259 overlaps B but not the primary. Locus identity
        // is anchored on the primary's span, so C stays a distinct locus rather than being chained in via B.
        SAMRecord record = newRecord(CHR_1, 1000, "100M");
        record.setMappingQuality(0);
        record.setAttribute("XA", CHR_1 + ",+1080,100M,0;" + CHR_1 + ",+1160,100M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.numLoci());
        assertEquals(0, result.updatedMapQuality());
    }

    @Test
    public void testRefOnlyMulti()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.numLoci());
    }

    @Test
    public void testTxPrimaryMulti()
    {
        ContigEntry entryB = new ContigEntry(
                "ensG_OTHER_T", 1, 100, "GO", "OTHER", "TO", "chr5", 1,
                List.of(new BaseRegion(2000, 2099)));
        List<ContigEntry> twoContigs = List.of(threeExonContig(), entryB);

        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "ensG_OTHER_T,+1,50M,0;");

        PlacementSelector selector = new PlacementSelector(twoContigs);
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.numLoci());
    }

    @Test
    public void testSupplementary()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "60M40H");
        record.setSupplementaryAlignmentFlag(true);

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertTrue(result.hasPlacement());
    }

    @Test
    public void testSupplementaryOnTxContigGetsLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setSupplementaryAlignmentFlag(true);
        record.setMappingQuality(0);

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertTrue(result.hasPlacement());
        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
        assertEquals(0, result.updatedMapQuality());
    }

    @Test
    public void testSupplementaryXaAlignmentsAreLiftedAndDeduplicatedAgainstSelf()
    {
        SAMRecord record = newRecord(CHR_1, 300, "50S50M");
        record.setSupplementaryAlignmentFlag(true);
        record.setAttribute(
                "XA", CHR_1 + ",+300,50S50M,0;" + TX_CONTIG + ",+151,50S50M,0;");

        LiftedRecord result = new PlacementSelector(contigMap()).resolve(record);

        assertEquals(2, result.liftedAlignments().size());
        assertEquals(CHR_1, result.liftedAlignments().get(1).LiftedChromosome);
        assertEquals(350, result.liftedAlignments().get(1).LiftedPos);
        assertEquals("50S50M", result.liftedAlignments().get(1).LiftedCigar);
        assertTrue(result.liftedAlignments().get(1).FromTxContig);
    }

    @Test
    public void testSupplementaryOnTxContigUnliftablePastEnd()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M"); // pos 251 past contigEnd(250)
        record.setSupplementaryAlignmentFlag(true);

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertFalse(result.hasPlacement());
        assertTrue(result.notes().contains("supp_translate_failed"));
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testPrimaryUnliftablePastEnd()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertFalse(result.hasPlacement());
        assertTrue(result.notes().contains("primary_translate_failed"));
    }

    @Test
    public void testPrimaryTrailingOverhangClampedToSoftClip()
    {
        SAMRecord record = newRecord(TX_CONTIG, 200, "100M"); // 49bp past contigEnd(250) -> trailing 49S

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertTrue(result.finalCigar().endsWith("49S"));
    }

    @Test
    public void testIntronRetRefBetterLeadingSoftClipBoundary()
    {
        // Tx alt at contig pos 101 carries its leading 20S at the exon-1/exon-2 boundary.
        SAMRecord record = newRecord(CHR_1, 300, "30M");
        record.setAttribute("XA", TX_CONTIG + ",+101,20S30M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testXaDedupDropsDuplicateAlts()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // two identical XA entries -> one alt retained
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,+5000,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + one deduped alt
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaDedupKeepsAltMatchingSelfButDropsXaDuplicate()
    {
        // Two XA entries lift to the same (chr1,100,50M): one Tx, one ref. XA dedup is XA-internal only, so the Tx alt
        // is kept and drives CONCORDANT while the duplicate ref XA collapses.
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;" + CHR_1 + ",+100,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + Tx alt; ref XA collapsed
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaWithMalformedNmStillLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // a garbled NM field must not silently drop the alt
        record.setAttribute("XA", "chr5,+5000,50M,not_a_number;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.liftedAlignments().size());
    }

    @Test
    public void testCrossLocusBothSplicedRemainsMultiLocus()
    {
        // Tx primary and ref alt are both spliced but at different loci: ambiguous, so the primary must not swap.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", "chr5,+5000,50M100N50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertTrue(result.primaryAlignment().FromTxContig);
    }

    // numLoci must reflect the deduped genomic-locus count, not the XA entry count (NH is derived from it).
    @Test
    public void testNumLociDedupesIdenticalLiftedXaEntries()
    {
        SAMRecord primary = newRecord(TX_CONTIG, 51, "100M"); // lifts to chr1:150 50M100N50M
        // four XA entries all lifting to the same locus -> numLoci still 1
        primary.setAttribute(
                "XA",
                TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(primary);

        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testNumLociCountsDistinctLiftedLoci()
    {
        SAMRecord primary = newRecord(CHR_1, 1000, "150M");
        primary.setAttribute(
                "XA",
                CHR_1 + ",+2000,150M,0;"
                        + CHR_1 + ",+3000,150M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(primary);

        assertEquals(3, result.numLoci());
    }

    // Hidden tie (XS==AS) on a ref-only primary outside any indexed exon: no evidence either way, so MAPQ holds at 0
    // because the equal-scoring alt bwa did not emit may be real.
    @Test
    public void testHiddenTieRefOnlyOutsideIndexedExonHoldsAtZero() throws Exception
    {
        SAMRecord record = newRecord(CHR_1, 1500, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        PlacementSelector noIndex = new PlacementSelector(contigMap());
        assertEquals(0, noIndex.resolve(record).updatedMapQuality());

        EnsemblAnnotationIndex annotationIndex = exonRegionIndex(CHR_1, List.of(new int[] { 1400, 1700 }));
        PlacementSelector withIndex = new PlacementSelector(contigMap(), annotationIndex);
        LiftedRecord result = withIndex.resolve(record);
        assertEquals(0, result.updatedMapQuality());
    }

    // Hidden tie with the primary outside any annotated exon: rescue stays blocked.
    @Test
    public void testHiddenTieOutsideExonKeepsMapQualityZero() throws Exception
    {
        SAMRecord record = newRecord(CHR_1, 5000, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        EnsemblAnnotationIndex annotationIndex = exonRegionIndex(
                CHR_1, List.of(new int[] { 1400, 1700 })); // primary at 5000 is intergenic
        PlacementSelector selector = new PlacementSelector(contigMap(), annotationIndex);
        assertEquals(0, selector.resolve(record).updatedMapQuality());
    }

    // decidePrimaryMapQuality args: (inputMapQuality, numLoci, hiddenTie, primaryFromTxContig, primaryInAnnotatedExon, randomTie).
    @Test
    public void testMapQualityPolicy_singleLocusZeroRescues()
    {
        assertEquals(60, PlacementSelector.decidePrimaryMapQuality(0, 1, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_hiddenTieRefPrimaryNoExonHoldsAtZero()
    {
        assertEquals(0, PlacementSelector.decidePrimaryMapQuality(0, 1, true, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_hiddenTieTxPrimaryRescues()
    {
        assertEquals(60, PlacementSelector.decidePrimaryMapQuality(0, 1, true, true, false, false));
    }

    @Test
    public void testMapQualityPolicy_hiddenTieInAnnotatedExonRescues()
    {
        assertEquals(60, PlacementSelector.decidePrimaryMapQuality(0, 1, true, false, true, false));
    }

    @Test
    public void testMapQualityPolicy_inputSixtyPassesAsRescued()
    {
        assertEquals(60, PlacementSelector.decidePrimaryMapQuality(60, 1, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_gradedMapQualityPassesThrough()
    {
        // a graded MAPQ is a real bwa signal, so it is left alone
        assertEquals(37, PlacementSelector.decidePrimaryMapQuality(37, 1, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_multiLocusNeverBumps()
    {
        assertEquals(0, PlacementSelector.decidePrimaryMapQuality(0, 2, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_randomTieNotBumped()
    {
        // a random-tie pick is a coin-flip among distinct placements, not a confident unique call
        assertEquals(0, PlacementSelector.decidePrimaryMapQuality(0, 1, false, false, false, true));
        assertEquals(60, PlacementSelector.decidePrimaryMapQuality(0, 1, false, false, false, false));
    }

    private static LiftedAlignment liftedAt(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, true, 0);
    }

    @Test
    public void testCountDistinctLociFromListRecountsPostExtension()
    {
        // The record overload backs the emit-time NH recompute: primary from primaryIndex, Dropped alts excluded, alts
        // overlapping the primary collapsed, so NH stays consistent with the XA tag.
        LiftedAlignment primary = liftedAt(CHR_1, 1000, "100M");   // span 1000-1099
        LiftedAlignment overlapping = liftedAt(CHR_1, 1050, "100M");   // 1050-1149 overlaps the primary
        LiftedAlignment distant = liftedAt(CHR_1, 5000, "100M");   // a genuinely distinct locus
        LiftedAlignment droppedDistant = liftedAt(CHR_1, 8000, "100M");
        droppedDistant.Dropped = true;

        assertEquals(1, countDistinctLociOf(primary));
        assertEquals(1, countDistinctLociOf(primary, overlapping));
        assertEquals(2, countDistinctLociOf(primary, distant));
        assertEquals(1, countDistinctLociOf(primary, droppedDistant));
        // no primary placement (an unmapped record's empty list) -> 1
        assertEquals(1, PlacementSelector.countDistinctLoci(LiftedRecord.unmapped("")));
    }

    private static int countDistinctLociOf(final LiftedAlignment... alignments)
    {
        return PlacementSelector.countDistinctLoci(recordBuilder().alignments(List.of(alignments)).build());
    }

    @Test
    public void testOppositeStrandXaAltsNotCollapsed()
    {
        // Two XA alts at the same locus and cigar but opposite strands are distinct placements, so the lifted dedup key
        // includes strand and both survive into the XA output.
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,-5000,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        String xa = result.xaTag();
        assertNotNull(xa);
        assertTrue("both strands kept: " + xa, xa.contains("chr5,+5000") && xa.contains("chr5,-5000"));
    }

    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";

    // Post-Step-1 the selector treats any surviving N as a real junction.
    private static final String TX_JUNCTION_CIGAR = "50M100N50M";
    private static final String FULL_MATCH_CIGAR = "100M";
    private static final String SOFTCLIP_CIGAR = "50M51S";

    private static LiftedAlignment tx(
            final String chrom, final int pos, final String cigar, final int numMismatches)
    {
        return new LiftedAlignment(chrom, pos, cigar, numMismatches, true, true, 1);
    }

    private static LiftedAlignment ref(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, true, 0);
    }

    private static LiftedAlignment refReverse(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, false, 0);
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

    private static boolean concordantOf(final LiftedAlignment... alignments)
    {
        return PlacementSelector.isConcordant(set(alignments));
    }

    // Concordant is the one evidence flag that changes the pick: it short-circuits apply() to keep bwa's primary.
    @Test
    public void testConcordant()
    {
        // one locus, ref and tx on the same gapless cigar: the two views agree
        assertTrue(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, 0)));

        // a surviving N means the two views disagree about splicing
        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, TX_JUNCTION_CIGAR, 0)));

        // same locus but different cigars
        assertFalse(concordantOf(
                ref(CHR1, 100, SOFTCLIP_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, 0)));

        // a single source cannot agree with itself, at one locus or several
        assertFalse(concordantOf(ref(CHR1, 100, FULL_MATCH_CIGAR)));
        assertFalse(concordantOf(tx(CHR1, 100, TX_JUNCTION_CIGAR, 0)));
        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                ref(CHR2, 200, FULL_MATCH_CIGAR)));

        // more than one locus is a contest even when the cigars match
        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR2, 200, FULL_MATCH_CIGAR, 0)));
    }

    @Test
    public void testConcordantSkipsGateDroppedAlt()
    {
        // The overhang gate marks a collapsed XA alt Dropped before the selector runs; ignoring it leaves a lone
        // ref source, so the otherwise-agreeing pair is not concordant.
        LiftedAlignment droppedTx = tx(CHR1, 100, FULL_MATCH_CIGAR, 0);
        droppedTx.Dropped = true;

        assertFalse(concordantOf(ref(CHR1, 100, FULL_MATCH_CIGAR), droppedTx));
    }

    // One locus, a soft-clipped ref and a contiguous tx placement that only score can separate.
    private static List<LiftedAlignment> contestedSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "50M51S");
        LiftedAlignment txAlt = tx(CHR1, 100, "100M", 0);
        return set(self, txAlt);
    }

    // Two loci: a ref and a tx placement.
    private static List<LiftedAlignment> multiLocusSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "151M");
        LiftedAlignment altB = tx(CHR2, 200, "151M", 0);
        return set(self, altB);
    }

    @Test
    public void testBwaPriorityKeepsSelf()
    {
        // MAPQ > 0: bwa ranked the placements, so its primary is kept regardless of score.
        List<LiftedAlignment> alignments = contestedSet();
        alignments.get(0).GenomicScore = 10;
        alignments.get(1).GenomicScore = 99;
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, alignments.get(0), 0, true);
        assertSame(alignments.get(0), outcome.effectivePrimary());
        assertEquals("", outcome.note());
    }

    @Test
    public void testUnscoredKeepsSelf()
    {
        // No placement scored (a split read left for Step 3): keep bwa's primary and drop nothing.
        List<LiftedAlignment> alignments = contestedSet();
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, alignments.get(0), 0, false);
        assertSame(alignments.get(0), outcome.effectivePrimary());
        assertEquals("", outcome.note());
        assertFalse(alignments.get(0).Dropped);
        assertFalse(alignments.get(1).Dropped);
    }

    @Test
    public void testDecisiveScoreWinsRegardlessOfSeed()
    {
        // The higher genomic score wins outright, whatever the seed.
        List<LiftedAlignment> refWins = contestedSet();
        refWins.get(0).GenomicScore = 90;
        refWins.get(1).GenomicScore = 50;
        PlacementSelector.ApplyResult refResult = PlacementSelector.apply(refWins, false, refWins.get(0), 1, false);
        assertSame(refWins.get(0), refResult.effectivePrimary());
        assertEquals("score", refResult.note());

        List<LiftedAlignment> txWins = contestedSet();
        txWins.get(0).GenomicScore = 40;
        txWins.get(1).GenomicScore = 88;
        PlacementSelector.ApplyResult txResult = PlacementSelector.apply(txWins, false, txWins.get(0), 0, false);
        assertSame(txWins.get(1), txResult.effectivePrimary());
        assertEquals("score", txResult.note());
        assertEquals(1, txResult.primaryIndex());
        assertFalse("loser rides in XA, not dropped", txWins.get(0).Dropped);
    }

    @Test
    public void testScoreTieFallsToSeededRandom()
    {
        // Equal scores fall to the seeded pick. contestedSet order is [self, tx]: seed 0 -> self, seed 1 -> tx.
        List<LiftedAlignment> even = contestedSet();
        even.get(0).GenomicScore = 70;
        even.get(1).GenomicScore = 70;
        PlacementSelector.ApplyResult evenResult = PlacementSelector.apply(even, false, even.get(0), 0, false);
        assertSame(even.get(0), evenResult.effectivePrimary());
        assertEquals("random", evenResult.note());

        List<LiftedAlignment> odd = contestedSet();
        odd.get(0).GenomicScore = 70;
        odd.get(1).GenomicScore = 70;
        PlacementSelector.ApplyResult oddResult = PlacementSelector.apply(odd, false, odd.get(0), 1, false);
        assertSame(odd.get(1), oddResult.effectivePrimary());
        assertEquals("random", oddResult.note());
        assertFalse("tie loser rides in XA, not dropped", odd.get(0).Dropped);
    }

    @Test
    public void testScoreTieJunctionBeatsSoftClipAtSameLocus()
    {
        // A spliced and a soft-clip placement at the same locus, tied on score: the junction wins outright, whatever
        // the seed. Seed 0 over set order [soft-clip, junction] would pick the soft-clip on a plain random tie.
        LiftedAlignment softClip = ref(CHR1, 100, SOFTCLIP_CIGAR);
        LiftedAlignment junction = tx(CHR1, 100, TX_JUNCTION_CIGAR, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, softClip, 0, false);
        assertSame(junction, outcome.effectivePrimary());
        assertEquals("junction", outcome.note());
        assertEquals(alignments.indexOf(junction), outcome.primaryIndex());
        assertFalse("soft-clip loser rides in XA, not dropped", softClip.Dropped);
    }

    @Test
    public void testScoreTieJunctionAtDifferentLocusStaysRandom()
    {
        // Junction and soft-clip sit at different loci, so the same-locus rule does not fire and the tie is random.
        LiftedAlignment softClip = ref(CHR1, 100, SOFTCLIP_CIGAR);
        LiftedAlignment junction = tx(CHR2, 200, TX_JUNCTION_CIGAR, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, softClip, 0, false);
        assertSame("seed 0 -> first placement", softClip, outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testMultiLocusDecisiveScorePicksBestLocus()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 60;
        alignments.get(1).GenomicScore = 130;
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, alignments.get(0), 0, false);
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
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, alignments.get(0), 1, false);
        assertSame("seed 1 -> second placement", alignments.get(1), outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testScoreTieCollapsesIdenticalPlacements()
    {
        // Self and a tx alt lift to the same locus and CIGAR, a third tx alt is a distinct spliced placement. The
        // identical pair collapses, so seed 1 ties over two placements and lands on the spliced one, not the duplicate.
        LiftedAlignment self = ref(CHR1, 100, "100M");
        LiftedAlignment txSame = tx(CHR1, 100, "100M", 0);
        LiftedAlignment txSpliced = tx(CHR1, 100, TX_JUNCTION_CIGAR, 0);
        List<LiftedAlignment> alignments = set(self, txSame, txSpliced);
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        alignments.get(2).GenomicScore = 100;
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, self, 1, false);
        assertSame(txSpliced, outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testMultiLocusTieMateProximityPicksMateLocus()
    {
        // Tied loci CHR1:100 and CHR2:200 with the mate on CHR2 near 200: mate proximity beats the seed, which would
        // otherwise pick CHR1:100.
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR2, 250, "100M", false);
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, alignments.get(0), 0, false, mate);
        assertSame(alignments.get(1), outcome.effectivePrimary());
        assertEquals("mate", outcome.note());
        assertEquals(1, outcome.primaryIndex());
        assertFalse("tie loser rides in XA, not dropped", alignments.get(0).Dropped);
    }

    @Test
    public void testMultiLocusTiePicksClosestMatePlacement()
    {
        // Both candidates are within the broad 1 Mb transcript-span limit, but only CHR1:100 forms a
        // plausible short fragment with the mate. Treating proximity as a yes/no threshold leaves this
        // to the read-name seed and can create a false 500 kb discordant fragment.
        LiftedAlignment close = ref(CHR1, 100, "100M");
        LiftedAlignment distant = tx(CHR1, 500_000, "100M", 0);
        close.GenomicScore = 100;
        distant.GenomicScore = 100;
        List<LiftedAlignment> alignments = set(close, distant);
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR1, 250, "100M", false);

        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(
                alignments, false, close, 1, false, mate);

        assertSame(close, outcome.effectivePrimary());
        assertEquals("mate", outcome.note());
    }

    @Test
    public void testClosestMatePlacementBeatsSupplementarySupport()
    {
        // A MAPQ-0 supplementary must not pull the primary hundreds of kilobases away when an equally
        // scoring placement forms a much tighter pair.
        LiftedAlignment close = ref(CHR1, 100, "100M");
        LiftedAlignment distant = tx(CHR1, 500_000, "50M50S", 0).withSupplementaryMerge(
                500_000, "100M", List.of(1), List.of(), 0, 0);
        close.GenomicScore = 100;
        distant.GenomicScore = 100;
        List<LiftedAlignment> alignments = set(close, distant);
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR1, 250, "100M", false);

        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(
                alignments, false, close, 1, false, mate);

        assertSame(close, outcome.effectivePrimary());
        assertEquals("mate", outcome.note());
    }

    @Test
    public void testPairProximityPrecedesCombinedAlignmentScore()
    {
        LiftedAlignment close = ref(CHR1, 100, "52M");
        LiftedAlignment distant = ref(CHR1, 500_000, "52M");
        close.GenomicScore = 50;
        distant.GenomicScore = 100;

        LiftedAlignment mate = ref(CHR1, 250, "52M");
        mate.GenomicScore = 100;

        PlacementSelector.PairApplyResult outcome = PlacementSelector.applyPair(
                set(close, distant), false, close, false,
                set(mate), false, mate, false, 0);

        assertSame(close, outcome.first().effectivePrimary());
        assertSame(mate, outcome.second().effectivePrimary());
        assertEquals("mate", outcome.first().note());
    }

    @Test
    public void testPairSelectionPrefersDeletionBeforeShorterDuplication()
    {
        LiftedAlignment deletion = ref(CHR1, 100, "100M");
        LiftedAlignment duplication = ref(CHR1, 1100, "100M");
        LiftedAlignment mate = refReverse(CHR1, 1000, "100M");

        PlacementSelector.PairApplyResult outcome = PlacementSelector.applyPair(
                set(deletion, duplication), false, deletion, false,
                set(mate), false, mate, true, 0);

        assertSame(deletion, outcome.first().effectivePrimary());
    }

    @Test
    public void testPairSvPreferenceStopsPastOneMegabase()
    {
        LiftedAlignment distantDeletion = ref(CHR1, 100, "100M");
        LiftedAlignment localInversion = refReverse(CHR1, 1_999_000, "100M");
        LiftedAlignment mate = refReverse(CHR1, 2_000_000, "100M");

        PlacementSelector.PairApplyResult outcome = PlacementSelector.applyPair(
                set(distantDeletion, localInversion), false, distantDeletion, false,
                set(mate), false, mate, true, 0);

        assertSame(localInversion, outcome.first().effectivePrimary());
    }

    @Test
    public void testPairSvPreferenceIncludesOneMegabaseBoundary()
    {
        LiftedAlignment boundaryDeletion = ref(CHR1, 100, "100M");
        LiftedAlignment localInversion = refReverse(CHR1, 1_000_000, "100M");
        LiftedAlignment mate = refReverse(CHR1, 1_000_199, "100M");

        PlacementSelector.PairApplyResult outcome = PlacementSelector.applyPair(
                set(localInversion, boundaryDeletion), false, localInversion, false,
                set(mate), false, mate, true, 0);

        assertSame(boundaryDeletion, outcome.first().effectivePrimary());
    }

    @Test
    public void testPositiveMapQualityPairKeepsBwaPlacements()
    {
        LiftedAlignment bwa = ref(CHR1, 1100, "100M");
        LiftedAlignment alternative = ref(CHR1, 100, "100M");
        LiftedAlignment mate = refReverse(CHR1, 1000, "100M");

        PlacementSelector.PairApplyResult outcome = PlacementSelector.applyPair(
                set(bwa, alternative), false, bwa, true,
                set(mate), false, mate, true, 0);

        assertSame(bwa, outcome.first().effectivePrimary());
    }

    @Test
    public void testPairSelectionIsIndependentOfMateOrder()
    {
        LiftedAlignment firstChr1 = ref(CHR1, 100, "52M");
        LiftedAlignment firstChr2 = ref(CHR2, 100, "52M");
        LiftedAlignment secondChr1 = ref(CHR1, 250, "52M");
        LiftedAlignment secondChr2 = ref(CHR2, 250, "52M");
        firstChr1.GenomicScore = 50;
        firstChr2.GenomicScore = 100;
        secondChr1.GenomicScore = 100;
        secondChr2.GenomicScore = 50;

        PlacementSelector.PairApplyResult forward = PlacementSelector.applyPair(
                set(firstChr1, firstChr2), false, firstChr1, false,
                set(secondChr1, secondChr2), false, secondChr1, false, 0);
        PlacementSelector.PairApplyResult reverse = PlacementSelector.applyPair(
                set(secondChr1, secondChr2), false, secondChr1, false,
                set(firstChr1, firstChr2), false, firstChr1, false, 0);

        assertSame(firstChr1, forward.first().effectivePrimary());
        assertSame(secondChr1, forward.second().effectivePrimary());
        assertSame(secondChr1, reverse.first().effectivePrimary());
        assertSame(firstChr1, reverse.second().effectivePrimary());
    }

    @Test
    public void testMateDistanceExcludesIntronicCigarSpan()
    {
        LiftedAlignment spliced = ref(CHR1, 100, "10M1000N10M");
        LiftedAlignment insideIntron = ref(CHR1, 500, "10M");
        LiftedAlignment overlappingExon = ref(CHR1, 1115, "10M");

        assertEquals(391, spliced.alignedBlockDistance(insideIntron));
        assertEquals(0, spliced.alignedBlockDistance(overlappingExon));
        assertEquals(Integer.MAX_VALUE, spliced.alignedBlockDistance(ref(CHR2, 100, "10M")));
    }

    @Test
    public void testMultiLocusTieMateOnNeitherChromStaysRandom()
    {
        // The mate is on a third chromosome, so proximity does not discriminate and the seed decides.
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt("chr9", 500, "100M", false);
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, alignments.get(0), 1, false, mate);
        assertSame("seed 1 -> second placement", alignments.get(1), outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testMultiLocusTieMateTooFarStaysRandom()
    {
        // The mate is on CHR1 but more than 1 Mb from the CHR1 locus, so it is not proximal.
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR1, 5_000_000, "100M", false);
        PlacementSelector.ApplyResult outcome = PlacementSelector.apply(alignments, false, alignments.get(0), 1, false, mate);
        assertSame("seed 1 -> second placement", alignments.get(1), outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }
}
