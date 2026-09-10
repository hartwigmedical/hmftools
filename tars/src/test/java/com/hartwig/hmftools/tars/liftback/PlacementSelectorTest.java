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

public class PlacementSelectorTest
{
    private static List<ContigEntry> contigMap()
    {
        return List.of(threeExonContig());
    }

    private static SAMRecord newRecord(final String contig, final int pos, final String cigar)
    {
        return TarsTestFixtures.mappedRecord("read", contig, pos, cigar);
    }

    private static SAMRecord newUnmappedRecord()
    {
        return TarsTestFixtures.unmappedRecord("read");
    }

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
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");

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
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");

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
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testIntronRetRefBetter()
    {
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,30M20S,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxSoftClipNotAtBoundaryFallsToAmbiguous()
    {
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,25M25S,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testMultiLocusTwoLoci()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.numLoci());
    }

    @Test
    public void testDistinctLocusAltBlocksMapQualityRescue()
    {
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
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setMappingQuality(0);
        record.setAttribute("XA", TX_CONTIG + ",+101,50S50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
        assertEquals(60, result.updatedMapQuality());
    }

    @Test
    public void testChainedOverlapAltNotMergedThroughPrimary()
    {
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
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");
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
        SAMRecord record = newRecord(TX_CONTIG, 200, "100M");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertTrue(result.finalCigar().endsWith("49S"));
    }

    @Test
    public void testIntronRetRefBetterLeadingSoftClipBoundary()
    {
        SAMRecord record = newRecord(CHR_1, 300, "30M");
        record.setAttribute("XA", TX_CONTIG + ",+101,20S30M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testXaDedupDropsDuplicateAlts()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,+5000,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.liftedAlignments().size());
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaDedupKeepsAltMatchingSelfButDropsXaDuplicate()
    {
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;" + CHR_1 + ",+100,50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.liftedAlignments().size());
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaWithMalformedNmStillLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,not_a_number;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertEquals(2, result.liftedAlignments().size());
    }

    @Test
    public void testCrossLocusBothSplicedRemainsMultiLocus()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", "chr5,+5000,50M100N50M,0;");

        PlacementSelector selector = new PlacementSelector(contigMap());
        LiftedRecord result = selector.resolve(record);

        assertTrue(result.primaryAlignment().FromTxContig);
    }

    @Test
    public void testNumLociDedupesIdenticalLiftedXaEntries()
    {
        SAMRecord primary = newRecord(TX_CONTIG, 51, "100M");
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

    @Test
    public void testHiddenTieOutsideExonKeepsMapQualityZero() throws Exception
    {
        SAMRecord record = newRecord(CHR_1, 5000, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        EnsemblAnnotationIndex annotationIndex = exonRegionIndex(
                CHR_1, List.of(new int[] { 1400, 1700 }));
        PlacementSelector selector = new PlacementSelector(contigMap(), annotationIndex);
        assertEquals(0, selector.resolve(record).updatedMapQuality());
    }

    @Test
    public void testMapQualityPolicy()
    {
        List<MapQualityCase> cases = List.of(
                new MapQualityCase("single locus", 0, 1, false, false, false, false, 60),
                new MapQualityCase("hidden tie", 0, 1, true, false, false, false, 0),
                new MapQualityCase("hidden tie on tx", 0, 1, true, true, false, false, 60),
                new MapQualityCase("hidden tie in exon", 0, 1, true, false, true, false, 60),
                new MapQualityCase("MAPQ 60", 60, 1, false, false, false, false, 60),
                new MapQualityCase("graded MAPQ", 37, 1, false, false, false, false, 37),
                new MapQualityCase("multiple loci", 0, 2, false, false, false, false, 0),
                new MapQualityCase("random tie", 0, 1, false, false, false, true, 0));

        for(MapQualityCase test : cases)
        {
            assertEquals(test.name(), test.expected(), PlacementSelector.decidePrimaryMapQuality(
                    test.input(), test.loci(), test.hiddenTie(), test.fromTx(), test.inExon(), test.randomTie()));
        }
    }

    private record MapQualityCase(
            String name, int input, int loci, boolean hiddenTie, boolean fromTx,
            boolean inExon, boolean randomTie, int expected)
    {
    }

    private static LiftedAlignment liftedAt(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, true, 0);
    }

    @Test
    public void testCountDistinctLociFromListRecountsPostExtension()
    {
        LiftedAlignment primary = liftedAt(CHR_1, 1000, "100M");
        LiftedAlignment overlapping = liftedAt(CHR_1, 1050, "100M");
        LiftedAlignment distant = liftedAt(CHR_1, 5000, "100M");
        LiftedAlignment droppedDistant = liftedAt(CHR_1, 8000, "100M");
        droppedDistant.Dropped = true;

        assertEquals(1, countDistinctLociOf(primary));
        assertEquals(1, countDistinctLociOf(primary, overlapping));
        assertEquals(2, countDistinctLociOf(primary, distant));
        assertEquals(1, countDistinctLociOf(primary, droppedDistant));
        assertEquals(1, PlacementSelector.countDistinctLoci(LiftedRecord.unmapped("")));
    }

    private static int countDistinctLociOf(final LiftedAlignment... alignments)
    {
        return PlacementSelector.countDistinctLoci(recordBuilder().alignments(List.of(alignments)).build());
    }

    @Test
    public void testOppositeStrandXaAltsNotCollapsed()
    {
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

    @Test
    public void testConcordant()
    {
        assertTrue(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, 0)));

        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, TX_JUNCTION_CIGAR, 0)));

        assertFalse(concordantOf(
                ref(CHR1, 100, SOFTCLIP_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, 0)));

        assertFalse(concordantOf(ref(CHR1, 100, FULL_MATCH_CIGAR)));
        assertFalse(concordantOf(tx(CHR1, 100, TX_JUNCTION_CIGAR, 0)));
        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                ref(CHR2, 200, FULL_MATCH_CIGAR)));

        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR2, 200, FULL_MATCH_CIGAR, 0)));
    }

    @Test
    public void testConcordantSkipsGateDroppedAlt()
    {
        LiftedAlignment droppedTx = tx(CHR1, 100, FULL_MATCH_CIGAR, 0);
        droppedTx.Dropped = true;

        assertFalse(concordantOf(ref(CHR1, 100, FULL_MATCH_CIGAR), droppedTx));
    }

    private static List<LiftedAlignment> contestedSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "50M51S");
        LiftedAlignment txAlt = tx(CHR1, 100, "100M", 0);
        return set(self, txAlt);
    }

    private static List<LiftedAlignment> multiLocusSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "151M");
        LiftedAlignment altB = tx(CHR2, 200, "151M", 0);
        return set(self, altB);
    }

    @Test
    public void testBwaPriorityKeepsSelf()
    {
        List<LiftedAlignment> alignments = contestedSet();
        alignments.get(0).GenomicScore = 10;
        alignments.get(1).GenomicScore = 99;
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, alignments.get(0), 0, true);
        assertSame(alignments.get(0), outcome.alignment());
        assertEquals("", outcome.reason());
    }

    @Test
    public void testUnscoredKeepsSelf()
    {
        List<LiftedAlignment> alignments = contestedSet();
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, alignments.get(0), 0, false);
        assertSame(alignments.get(0), outcome.alignment());
        assertEquals("", outcome.reason());
        assertFalse(alignments.get(0).Dropped);
        assertFalse(alignments.get(1).Dropped);
    }

    @Test
    public void testDecisiveScoreWinsRegardlessOfSeed()
    {
        List<LiftedAlignment> refWins = contestedSet();
        refWins.get(0).GenomicScore = 90;
        refWins.get(1).GenomicScore = 50;
        PlacementSelector.Selection refResult = PlacementSelector.select(refWins, false, refWins.get(0), 1, false);
        assertSame(refWins.get(0), refResult.alignment());
        assertEquals("score", refResult.reason());

        List<LiftedAlignment> txWins = contestedSet();
        txWins.get(0).GenomicScore = 40;
        txWins.get(1).GenomicScore = 88;
        PlacementSelector.Selection txResult = PlacementSelector.select(txWins, false, txWins.get(0), 0, false);
        assertSame(txWins.get(1), txResult.alignment());
        assertEquals("score", txResult.reason());
        assertEquals(1, txResult.alignmentIndex());
        assertFalse("loser rides in XA, not dropped", txWins.get(0).Dropped);
    }

    @Test
    public void testScoreTieFallsToSeededRandom()
    {
        List<LiftedAlignment> even = contestedSet();
        even.get(0).GenomicScore = 70;
        even.get(1).GenomicScore = 70;
        PlacementSelector.Selection evenResult = PlacementSelector.select(even, false, even.get(0), 0, false);
        assertSame(even.get(0), evenResult.alignment());
        assertEquals("random", evenResult.reason());

        List<LiftedAlignment> odd = contestedSet();
        odd.get(0).GenomicScore = 70;
        odd.get(1).GenomicScore = 70;
        PlacementSelector.Selection oddResult = PlacementSelector.select(odd, false, odd.get(0), 1, false);
        assertSame(odd.get(1), oddResult.alignment());
        assertEquals("random", oddResult.reason());
        assertFalse("tie loser rides in XA, not dropped", odd.get(0).Dropped);
    }

    @Test
    public void testScoreTieJunctionBeatsSoftClipAtSameLocus()
    {
        LiftedAlignment softClip = ref(CHR1, 100, SOFTCLIP_CIGAR);
        LiftedAlignment junction = tx(CHR1, 100, TX_JUNCTION_CIGAR, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, softClip, 0, false);
        assertSame(junction, outcome.alignment());
        assertEquals("junction", outcome.reason());
        assertEquals(alignments.indexOf(junction), outcome.alignmentIndex());
        assertFalse("soft-clip loser rides in XA, not dropped", softClip.Dropped);
    }

    @Test
    public void testScoreTieJunctionAtDifferentLocusStaysRandom()
    {
        LiftedAlignment softClip = ref(CHR1, 100, SOFTCLIP_CIGAR);
        LiftedAlignment junction = tx(CHR2, 200, TX_JUNCTION_CIGAR, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, softClip, 0, false);
        assertSame("seed 0 -> first placement", softClip, outcome.alignment());
        assertEquals("random", outcome.reason());
    }

    @Test
    public void testMultiLocusDecisiveScorePicksBestLocus()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 60;
        alignments.get(1).GenomicScore = 130;
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, alignments.get(0), 0, false);
        assertSame(alignments.get(1), outcome.alignment());
        assertEquals("score", outcome.reason());
        assertFalse("all placements ride in XA", alignments.get(0).Dropped);
    }

    @Test
    public void testMultiLocusScoreTieSeededRandom()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, alignments.get(0), 1, false);
        assertSame("seed 1 -> second placement", alignments.get(1), outcome.alignment());
        assertEquals("random", outcome.reason());
    }

    @Test
    public void testScoreTieCollapsesIdenticalPlacements()
    {
        LiftedAlignment self = ref(CHR1, 100, "100M");
        LiftedAlignment txSame = tx(CHR1, 100, "100M", 0);
        LiftedAlignment txSpliced = tx(CHR1, 100, TX_JUNCTION_CIGAR, 0);
        List<LiftedAlignment> alignments = set(self, txSame, txSpliced);
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        alignments.get(2).GenomicScore = 100;
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, self, 1, false);
        assertSame(txSpliced, outcome.alignment());
        assertEquals("random", outcome.reason());
    }

    @Test
    public void testMultiLocusTieMateProximityPicksMateLocus()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR2, 250, "100M", false);
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, alignments.get(0), 0, false, mate);
        assertSame(alignments.get(1), outcome.alignment());
        assertEquals("mate", outcome.reason());
        assertEquals(1, outcome.alignmentIndex());
        assertFalse("tie loser rides in XA, not dropped", alignments.get(0).Dropped);
    }

    @Test
    public void testMultiLocusTiePicksClosestMatePlacement()
    {
        LiftedAlignment close = ref(CHR1, 100, "100M");
        LiftedAlignment distant = tx(CHR1, 500_000, "100M", 0);
        close.GenomicScore = 100;
        distant.GenomicScore = 100;
        List<LiftedAlignment> alignments = set(close, distant);
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR1, 250, "100M", false);

        PlacementSelector.Selection outcome = PlacementSelector.select(
                alignments, false, close, 1, false, mate);

        assertSame(close, outcome.alignment());
        assertEquals("mate", outcome.reason());
    }

    @Test
    public void testClosestMatePlacementBeatsSupplementarySupport()
    {
        LiftedAlignment close = ref(CHR1, 100, "100M");
        LiftedAlignment distant = tx(CHR1, 500_000, "50M50S", 0).withSupplementaryMerge(
                500_000, "100M", List.of(1), List.of(), 0, 0);
        close.GenomicScore = 100;
        distant.GenomicScore = 100;
        List<LiftedAlignment> alignments = set(close, distant);
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR1, 250, "100M", false);

        PlacementSelector.Selection outcome = PlacementSelector.select(
                alignments, false, close, 1, false, mate);

        assertSame(close, outcome.alignment());
        assertEquals("mate", outcome.reason());
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

        PlacementSelector.PairSelection outcome = PlacementSelector.selectPair(
                set(close, distant), false, close, false,
                set(mate), false, mate, false, 0);

        assertSame(close, outcome.first().alignment());
        assertSame(mate, outcome.second().alignment());
        assertEquals("mate", outcome.first().reason());
    }

    @Test
    public void testPairSelectionPrefersDeletionBeforeShorterDuplication()
    {
        LiftedAlignment deletion = ref(CHR1, 100, "100M");
        LiftedAlignment duplication = ref(CHR1, 1100, "100M");
        LiftedAlignment mate = refReverse(CHR1, 1000, "100M");

        PlacementSelector.PairSelection outcome = PlacementSelector.selectPair(
                set(deletion, duplication), false, deletion, false,
                set(mate), false, mate, true, 0);

        assertSame(deletion, outcome.first().alignment());
    }

    @Test
    public void testPairSvPreferenceStopsPastOneMegabase()
    {
        LiftedAlignment distantDeletion = ref(CHR1, 100, "100M");
        LiftedAlignment localInversion = refReverse(CHR1, 1_999_000, "100M");
        LiftedAlignment mate = refReverse(CHR1, 2_000_000, "100M");

        PlacementSelector.PairSelection outcome = PlacementSelector.selectPair(
                set(distantDeletion, localInversion), false, distantDeletion, false,
                set(mate), false, mate, true, 0);

        assertSame(localInversion, outcome.first().alignment());
    }

    @Test
    public void testPairSvPreferenceIncludesOneMegabaseBoundary()
    {
        LiftedAlignment boundaryDeletion = ref(CHR1, 100, "100M");
        LiftedAlignment localInversion = refReverse(CHR1, 1_000_000, "100M");
        LiftedAlignment mate = refReverse(CHR1, 1_000_199, "100M");

        PlacementSelector.PairSelection outcome = PlacementSelector.selectPair(
                set(localInversion, boundaryDeletion), false, localInversion, false,
                set(mate), false, mate, true, 0);

        assertSame(boundaryDeletion, outcome.first().alignment());
    }

    @Test
    public void testPositiveMapQualityPairKeepsBwaPlacements()
    {
        LiftedAlignment bwa = ref(CHR1, 1100, "100M");
        LiftedAlignment alternative = ref(CHR1, 100, "100M");
        LiftedAlignment mate = refReverse(CHR1, 1000, "100M");

        PlacementSelector.PairSelection outcome = PlacementSelector.selectPair(
                set(bwa, alternative), false, bwa, true,
                set(mate), false, mate, true, 0);

        assertSame(bwa, outcome.first().alignment());
    }

    @Test
    public void testPairSelectionFallsBackWhenAllPlacementsAreDropped()
    {
        LiftedAlignment first = ref(CHR1, 100, "100M");
        LiftedAlignment second = refReverse(CHR1, 300, "100M");
        first.Dropped = true;
        second.Dropped = true;

        PlacementSelector.PairSelection outcome = PlacementSelector.selectPair(
                set(first), false, first, false,
                set(second), false, second, false, 0);

        assertSame(first, outcome.first().alignment());
        assertSame(second, outcome.second().alignment());
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

        PlacementSelector.PairSelection forward = PlacementSelector.selectPair(
                set(firstChr1, firstChr2), false, firstChr1, false,
                set(secondChr1, secondChr2), false, secondChr1, false, 0);
        PlacementSelector.PairSelection reverse = PlacementSelector.selectPair(
                set(secondChr1, secondChr2), false, secondChr1, false,
                set(firstChr1, firstChr2), false, firstChr1, false, 0);

        assertSame(firstChr1, forward.first().alignment());
        assertSame(secondChr1, forward.second().alignment());
        assertSame(secondChr1, reverse.first().alignment());
        assertSame(firstChr1, reverse.second().alignment());
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
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt("chr9", 500, "100M", false);
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, alignments.get(0), 1, false, mate);
        assertSame("seed 1 -> second placement", alignments.get(1), outcome.alignment());
        assertEquals("random", outcome.reason());
    }

    @Test
    public void testMultiLocusTieMateTooFarStaysRandom()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR1, 5_000_000, "100M", false);
        PlacementSelector.Selection outcome = PlacementSelector.select(alignments, false, alignments.get(0), 1, false, mate);
        assertSame("seed 1 -> second placement", alignments.get(1), outcome.alignment());
        assertEquals("random", outcome.reason());
    }
}
