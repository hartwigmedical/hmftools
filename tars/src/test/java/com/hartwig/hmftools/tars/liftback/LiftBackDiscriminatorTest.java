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

// Lift + decide over a record: coordinate lifting and result assembly, the primary pick over a hand-built
// alignment set, and the MAPQ policy. Cigar translation itself is covered in ContigTranslatorTest.
public class LiftBackDiscriminatorTest
{
    private static List<ContigEntry> contigMap()
    {
        return List.of(threeExonContig());
    }

    // no pairing flags: these tests must not carry them, and the read name seeds the random tie-break.
    private static SAMRecord newRecord(final String contig, final int pos, final String cigar)
    {
        return TarsTestFixtures.mappedRecord("read", contig, pos, cigar);
    }

    private static SAMRecord newUnmappedRecord()
    {
        return TarsTestFixtures.unmappedRecord("read");
    }

    // Ultima single-end: no pairing flags at all, so any unguarded htsjdk pair getter on the resolve path throws.
    // The lifted placement must match the paired case exactly - pairing plays no part in the decision.
    @Test
    public void testSingleEndPrimaryLifted()
    {
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(TarsTestFixtures.unpairedPrimaryRecord(TX_CONTIG, 51, "100M"));

        assertTrue(result.hasPlacement());
        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
    }

    @Test
    public void testSingleEndSupplementaryLifted()
    {
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(
                TarsTestFixtures.unpairedSupplementaryRecord(TX_CONTIG, 51, "100M", TX_CONTIG + ",1,+,60M40S,60,0;"));

        assertTrue(result.hasPlacement());
        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
    }

    @Test
    public void testUnmapped()
    {
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(newUnmappedRecord());

        assertFalse(result.hasPlacement());
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testRefOnlyExonic()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "150M");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

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

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(100, result.finalPos());
        assertEquals("50M", result.finalCigar());
        assertFalse(result.hasNCigar());
    }

    @Test
    public void testTxPrimaryUniqueJunctionCrosser()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M"); // crosses exon 1 -> exon 2 (intron 200-299)

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
    }

    @Test
    public void testRefTxAgree()
    {
        SAMRecord record = newRecord(CHR_1, 100, "50M"); // ref primary; Tx alt at contig pos 1 -> same locus+CIGAR
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testIntronRetRefBetter()
    {
        // ref chr1:170 50M (last 20bp in intron); Tx pos 71 30M20S lifts to chr1:170-199 soft-clipped at exon boundary.
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,30M20S,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxSoftClipNotAtBoundaryFallsToAmbiguous()
    {
        // ref 50M full match; Tx 25M25S - 25M ends mid-exon (chr1:194), trailing clip NOT at exon boundary -> AMBIGUOUS.
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,25M25S,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testMultiLocusTwoLoci()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // Tx primary + ref alt on different chrom -> two loci
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(2, result.numLoci());
    }

    @Test
    public void testDistinctLocusAltBlocksMapQualityRescue()
    {
        // Tx primary (perfect) + a lower-scoring XA alt (3 mismatches) at a distinct locus. bwa left the read
        // MAPQ 0; two distinct genomic loci means it stays a multimapper - TARS does not override bwa's call with a
        // weaker reconstructed score - so numLoci == 2 and MAPQ is held at 0.
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setMappingQuality(0);
        record.setAttribute("XA", "chr5,+5000,50M,3;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(2, result.numLoci());
        assertEquals(0, result.updatedMapQuality());
    }

    @Test
    public void testNestedSpanAltCountsAsSingleLocusAndBumpsMapQuality()
    {
        // Regression: a 5'-softclipped isoform alt begins at a downstream exon of the SAME placement. It lifts to a
        // different start (chr1:300) but a span (300-349) nested inside the junction-crossing primary's (chr1:150-349),
        // so it is the same genomic locus, not a second one. A MAPQ-0 read here is therefore a single-locus placement
        // and must bump to 60 -- keying loci on exact start (the old behaviour) miscounted this as two loci and held it at 0.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M"); // -> chr1:150 50M100N50M, genomic span 150-349
        record.setMappingQuality(0);
        record.setAttribute("XA", TX_CONTIG + ",+101,50S50M,0;"); // -> chr1:300 50S50M, span 300-349 (nested)

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(1, result.numLoci());
        assertEquals(60, result.updatedMapQuality());
    }

    @Test
    public void testChainedOverlapAltNotMergedThroughPrimary()
    {
        // A tandem-repeat-style read: primary chr1:1000-1099; alt B chr1:1080-1179 overlaps the primary; alt C
        // chr1:1160-1259 overlaps B but NOT the primary. C is a genuinely distinct locus, so numLoci must be 2 and a
        // MAPQ-0 read stays 0. Locus identity is anchored on the PRIMARY's span - C must not be chained back into the
        // primary via B (a naive interval-merge over all alignments would wrongly collapse to one locus and bump to 60).
        SAMRecord record = newRecord(CHR_1, 1000, "100M");
        record.setMappingQuality(0);
        record.setAttribute("XA", CHR_1 + ",+1080,100M,0;" + CHR_1 + ",+1160,100M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(2, result.numLoci());
        assertEquals(0, result.updatedMapQuality());
    }

    @Test
    public void testRefOnlyMulti()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "50M"); // ref primary + ref alt on different chrom, no tx
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

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

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(twoContigs);
        LiftedRecord result = resolver.resolve(record);

        assertEquals(2, result.numLoci());
    }

    @Test
    public void testSupplementary()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "60M40H");
        record.setSupplementaryAlignmentFlag(true);

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertTrue(result.hasPlacement());
    }

    @Test
    public void testSupplementaryOnTxContigGetsLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setSupplementaryAlignmentFlag(true);

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertTrue(result.hasPlacement());
        assertEquals(CHR_1, result.finalChromosome());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
    }

    @Test
    public void testSupplementaryOnTxContigUnliftablePastEnd()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M"); // pos 251 past contigEnd(250)
        record.setSupplementaryAlignmentFlag(true);

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertFalse(result.hasPlacement());
        assertTrue(result.notes().contains("supp_translate_failed"));
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testPrimaryUnliftablePastEnd()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertFalse(result.hasPlacement());
        assertTrue(result.notes().contains("primary_translate_failed"));
    }

    @Test
    public void testPrimaryTrailingOverhangClampedToSoftClip()
    {
        SAMRecord record = newRecord(TX_CONTIG, 200, "100M"); // 49bp past contigEnd(250) -> trailing 49S

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertTrue(result.finalCigar().endsWith("49S"));
    }

    @Test
    public void testIntronRetRefBetterLeadingSoftClipBoundary()
    {
        // Tx alt contig pos 101 with leading 20S at exon-1/exon-2 boundary; ref chr1:300 full match.
        SAMRecord record = newRecord(CHR_1, 300, "30M");
        record.setAttribute("XA", TX_CONTIG + ",+101,20S30M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(1, result.numLoci());
    }

    @Test
    public void testXaDedupDropsDuplicateAlts()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // two identical XA entries -> one alt retained
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,+5000,50M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + one deduped alt
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaDedupKeepsAltMatchingSelfButDropsXaDuplicate()
    {
        // Two XA entries lifting to the same (chr1,100,50M): one Tx, one ref. XA dedup is XA-internal only,
        // so the Tx alt is kept (drives CONCORDANT); the duplicate ref XA collapses.
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;" + CHR_1 + ",+100,50M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + Tx alt; ref XA collapsed
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaWithMalformedNmStillLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // garbled NM field must not silently drop the alt
        record.setAttribute("XA", "chr5,+5000,50M,not_a_number;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size());
    }

    @Test
    public void testCrossLocusBothSplicedRemainsMultiLocus()
    {
        // Both Tx primary and ref XA alt have spliced CIGARs at different loci - ambiguous, must NOT swap.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", "chr5,+5000,50M100N50M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        assertTrue(result.primaryAlignment().FromTxContig);
    }

    @Test
    public void testStatsEndToEndSmoke()
    {
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftBackStatistics stats = new LiftBackStatistics();

        SAMRecord r1 = newRecord(CHR_1, 1000, "150M"); // ref only
        stats.record(r1, resolver.resolve(r1));

        SAMRecord r2 = newRecord(TX_CONTIG, 51, "100M"); // tx only, junction-crosser
        stats.record(r2, resolver.resolve(r2));

        SAMRecord r3 = newRecord(CHR_1, 100, "50M"); // ref and tx agreeing, MAPQ=0
        r3.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");
        r3.setMappingQuality(0);
        stats.record(r3, resolver.resolve(r3));

        SAMRecord r4 = newUnmappedRecord(); // UNMAPPED
        stats.record(r4, resolver.resolve(r4));

        assertEquals(4, stats.total());
        assertEquals(3, stats.resolved());
        // the fourth arrived unmapped, so nothing was lost by lift-back
        assertEquals(0, stats.unmapped());
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

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(primary);

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

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(primary);

        assertEquals(3, result.numLoci());
    }

    // Hidden tie (XS==AS) on a ref-only primary landing outside any indexed exon: no vouching evidence, so the
    // unresolved hidden tie holds MAPQ at 0 (the equal-scoring alt bwa did not emit may be real).
    @Test
    public void testHiddenTieRefOnlyOutsideIndexedExonHoldsAtZero() throws Exception
    {
        SAMRecord record = newRecord(CHR_1, 1500, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        LiftBackDiscriminator noIndex = new LiftBackDiscriminator(contigMap());
        assertEquals(0, noIndex.resolve(record).updatedMapQuality());

        ExonRegionIndex exonIndex = exonRegionIndex(CHR_1, List.of(new int[] { 1400, 1700 }));
        LiftBackDiscriminator withIndex = new LiftBackDiscriminator(contigMap(), exonIndex);
        LiftedRecord result = withIndex.resolve(record);
        assertEquals(0, result.updatedMapQuality());
    }

    // Hidden tie with primary outside any annotated exon: rescue stays blocked.
    @Test
    public void testHiddenTieOutsideExonKeepsMapQualityZero() throws Exception
    {
        SAMRecord record = newRecord(CHR_1, 5000, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        ExonRegionIndex exonIndex = exonRegionIndex(
                CHR_1, List.of(new int[] { 1400, 1700 })); // exon at 1400-1700; primary at 5000 is intergenic
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap(), exonIndex);
        assertEquals(0, resolver.resolve(record).updatedMapQuality());
    }

    // Direct unit tests for the extracted MAPQ policy; independent of LiftBackDiscriminator / SAMRecord plumbing.
    // decidePrimaryMapQuality positional args: (inputMapQuality, numLoci, hiddenTie, primaryFromTxContig, primaryInAnnotatedExon).
    @Test
    public void testMapQualityPolicy_singleLocusZeroRescues()
    {
        // single locus, MAPQ0, no hidden tie -> 60
        assertEquals(60, LiftBackDiscriminator.decidePrimaryMapQuality(0, 1, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_hiddenTieRefPrimaryNoExonHoldsAtZero()
    {
        // unresolved hidden tie holds at 0
        assertEquals(0, LiftBackDiscriminator.decidePrimaryMapQuality(0, 1, true, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_hiddenTieTxPrimaryRescues()
    {
        // tx provenance overrides hidden tie
        assertEquals(60, LiftBackDiscriminator.decidePrimaryMapQuality(0, 1, true, true, false, false));
    }

    @Test
    public void testMapQualityPolicy_hiddenTieInAnnotatedExonRescues()
    {
        // annotated exon overrides hidden tie
        assertEquals(60, LiftBackDiscriminator.decidePrimaryMapQuality(0, 1, true, false, true, false));
    }

    @Test
    public void testMapQualityPolicy_inputSixtyPassesAsRescued()
    {
        // input 60 passes through unchanged
        assertEquals(60, LiftBackDiscriminator.decidePrimaryMapQuality(60, 1, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_gradedMapQualityPassesThrough()
    {
        // graded signal; leave alone
        assertEquals(37, LiftBackDiscriminator.decidePrimaryMapQuality(37, 1, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_multiLocusNeverBumps()
    {
        // multi-locus never bumped
        assertEquals(0, LiftBackDiscriminator.decidePrimaryMapQuality(0, 2, false, false, false, false));
    }

    @Test
    public void testMapQualityPolicy_randomTieNotBumped()
    {
        // a random-tie pick is a coin-flip among distinct placements, so it is not a confident unique call
        assertEquals(0, LiftBackDiscriminator.decidePrimaryMapQuality(0, 1, false, false, false, true));
        assertEquals(60, LiftBackDiscriminator.decidePrimaryMapQuality(0, 1, false, false, false, false));
    }

    private static LiftedAlignment liftedAt(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, false, true, 0);
    }

    @Test
    public void testCountDistinctLociFromListRecountsPostExtension()
    {
        // The record overload backs the emit-time NH recompute: it takes the primary from primaryIndex, drops Dropped
        // alts, and collapses alts overlapping the primary - so NH stays consistent with the XA after
        // reconcileChosenPrimary's alt extension mutates the shared list.
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
        assertEquals(1, LiftBackDiscriminator.countDistinctLoci(LiftedRecord.unmapped("")));
    }

    private static int countDistinctLociOf(final LiftedAlignment... alignments)
    {
        return LiftBackDiscriminator.countDistinctLoci(recordBuilder().alignments(List.of(alignments)).build());
    }

    @Test
    public void testOppositeStrandXaAltsNotCollapsed()
    {
        // Two XA alts at the same genomic locus and cigar but opposite strands are distinct placements. The lifted
        // dedup key includes strand, so both survive into the XA output; previously they collapsed to one, silently
        // dropping a strand-distinct placement.
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,-5000,50M,0;");

        LiftBackDiscriminator resolver = new LiftBackDiscriminator(contigMap());
        LiftedRecord result = resolver.resolve(record);

        String xa = result.xaTag();
        assertNotNull(xa);
        assertTrue("both strands kept: " + xa, xa.contains("chr5,+5000") && xa.contains("chr5,-5000"));
    }

    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";

    // carries an N: post-Step-1 the discriminator treats any surviving N as a real junction.
    private static final String TX_JUNCTION_CIGAR = "50M100N50M";
    private static final String FULL_MATCH_CIGAR = "100M";
    private static final String SOFTCLIP_CIGAR = "50M51S";

    private static LiftedAlignment tx(
            final String chrom, final int pos, final String cigar, final boolean softClipAtBoundary, final int numMismatches)
    {
        return new LiftedAlignment(chrom, pos, cigar, numMismatches, true, softClipAtBoundary, true, 1);
    }

    private static LiftedAlignment ref(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, false, true, 0);
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
        return LiftBackDiscriminator.isConcordant(set(alignments));
    }

    // Concordant is the one evidence flag that changes the pick: it short-circuits apply() to keep bwa's primary.
    @Test
    public void testConcordant()
    {
        // one locus, ref and tx on the same gapless cigar: the two views agree, nothing to choose between
        assertTrue(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0)));

        // a surviving N means the two views disagree about splicing
        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0)));

        // same locus but different cigars
        assertFalse(concordantOf(
                ref(CHR1, 100, SOFTCLIP_CIGAR),
                tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0)));

        // a single source cannot agree with itself, at one locus or several
        assertFalse(concordantOf(ref(CHR1, 100, FULL_MATCH_CIGAR)));
        assertFalse(concordantOf(tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0)));
        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                ref(CHR2, 200, FULL_MATCH_CIGAR)));

        // more than one locus is a contest even when the cigars match
        assertFalse(concordantOf(
                ref(CHR1, 100, FULL_MATCH_CIGAR),
                tx(CHR2, 200, FULL_MATCH_CIGAR, false, 0)));
    }

    @Test
    public void testConcordantSkipsGateDroppedAlt()
    {
        // an XA alt the overhang gate collapsed to a contiguous alignment is marked Dropped before the discriminator
        // runs; ignoring it leaves a lone ref source, so the otherwise-agreeing pair is not concordant.
        LiftedAlignment droppedTx = tx(CHR1, 100, FULL_MATCH_CIGAR, false, 0);
        droppedTx.Dropped = true;

        assertFalse(concordantOf(ref(CHR1, 100, FULL_MATCH_CIGAR), droppedTx));
    }

    // ---- apply(): score-based primary pick (no shape rules) ----

    // single locus: a ref (softclipped) and a tx (contiguous) candidate that only score can separate.
    private static List<LiftedAlignment> contestedSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "50M51S");
        LiftedAlignment txAlt = tx(CHR1, 100, "100M", false, 0);
        return set(self, txAlt);
    }

    // two loci: a ref and a tx placement.
    private static List<LiftedAlignment> multiLocusSet()
    {
        LiftedAlignment self = ref(CHR1, 100, "151M");
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
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, alignments.get(0), 0, true);
        assertSame(alignments.get(0), outcome.effectivePrimary());
        assertEquals("", outcome.note());
    }

    @Test
    public void testUnscoredKeepsSelf()
    {
        // contested but no candidate scored (a split read left for Step 3): keep bwa's primary, drop nothing.
        List<LiftedAlignment> alignments = contestedSet();
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, alignments.get(0), 0, false);
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
        LiftBackDiscriminator.ApplyResult refResult = LiftBackDiscriminator.apply(refWins, false, refWins.get(0), 1, false);
        assertSame(refWins.get(0), refResult.effectivePrimary());
        assertEquals("score", refResult.note());

        List<LiftedAlignment> txWins = contestedSet();
        txWins.get(0).GenomicScore = 40;
        txWins.get(1).GenomicScore = 88;
        LiftBackDiscriminator.ApplyResult txResult = LiftBackDiscriminator.apply(txWins, false, txWins.get(0), 0, false);
        assertSame(txWins.get(1), txResult.effectivePrimary());
        assertEquals("score", txResult.note());
        assertEquals(1, txResult.primaryIndex());
        assertFalse("loser rides in XA, not dropped", txWins.get(0).Dropped);
    }

    @Test
    public void testScoreTieFallsToSeededRandom()
    {
        // equal scores -> seeded pick, reproducible. contestedSet order is [self, tx]; seed 0 -> self, seed 1 -> tx.
        List<LiftedAlignment> even = contestedSet();
        even.get(0).GenomicScore = 70;
        even.get(1).GenomicScore = 70;
        LiftBackDiscriminator.ApplyResult evenResult = LiftBackDiscriminator.apply(even, false, even.get(0), 0, false);
        assertSame(even.get(0), evenResult.effectivePrimary());
        assertEquals("random", evenResult.note());

        List<LiftedAlignment> odd = contestedSet();
        odd.get(0).GenomicScore = 70;
        odd.get(1).GenomicScore = 70;
        LiftBackDiscriminator.ApplyResult oddResult = LiftBackDiscriminator.apply(odd, false, odd.get(0), 1, false);
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
        LiftedAlignment junction = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, softClip, 0, false);
        assertSame(junction, outcome.effectivePrimary());
        assertEquals("junction", outcome.note());
        assertEquals(alignments.indexOf(junction), outcome.primaryIndex());
        assertFalse("soft-clip loser rides in XA, not dropped", softClip.Dropped);
    }

    @Test
    public void testScoreTieJunctionAtDifferentLocusStaysRandom()
    {
        // junction and soft-clip sit at different loci, so the same-locus rule does not fire and the tie is random.
        LiftedAlignment softClip = ref(CHR1, 100, SOFTCLIP_CIGAR);
        LiftedAlignment junction = tx(CHR2, 200, TX_JUNCTION_CIGAR, false, 0);
        List<LiftedAlignment> alignments = set(softClip, junction);
        alignments.get(0).GenomicScore = 80;
        alignments.get(1).GenomicScore = 80;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, softClip, 0, false);
        assertSame("seed 0 -> first candidate", softClip, outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }

    @Test
    public void testMultiLocusDecisiveScorePicksBestLocus()
    {
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 60;
        alignments.get(1).GenomicScore = 130;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, alignments.get(0), 0, false);
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
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, alignments.get(0), 1, false);
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
        LiftedAlignment txSame = tx(CHR1, 100, "100M", false, 0);
        LiftedAlignment txSpliced = tx(CHR1, 100, TX_JUNCTION_CIGAR, false, 0);
        List<LiftedAlignment> alignments = set(self, txSame, txSpliced);
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        alignments.get(2).GenomicScore = 100;
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, self, 1, false);
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
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR2, 250, "100M", false);
        // seed 0 would pick CHR1:100; mate proximity overrides to the CHR2 locus.
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, alignments.get(0), 0, false, mate);
        assertSame(alignments.get(1), outcome.effectivePrimary());
        assertEquals("mate", outcome.note());
        assertEquals(1, outcome.primaryIndex());
        assertFalse("tie loser rides in XA, not dropped", alignments.get(0).Dropped);
    }

    @Test
    public void testMultiLocusTieMateOnNeitherChromStaysRandom()
    {
        // the mate is on a third chromosome, so proximity does not discriminate and the seed decides.
        List<LiftedAlignment> alignments = multiLocusSet();
        alignments.get(0).GenomicScore = 100;
        alignments.get(1).GenomicScore = 100;
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt("chr9", 500, "100M", false);
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, alignments.get(0), 1, false, mate);
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
        LiftedRecord mate = TarsTestFixtures.liftedRecordAt(CHR1, 5_000_000, "100M", false);
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(alignments, false, alignments.get(0), 1, false, mate);
        assertSame("seed 1 -> second candidate", alignments.get(1), outcome.effectivePrimary());
        assertEquals("random", outcome.note());
    }
}
