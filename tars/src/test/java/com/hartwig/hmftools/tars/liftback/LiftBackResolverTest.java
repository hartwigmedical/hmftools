package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.exonRegionIndex;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// Per-category tests for LiftBackResolver plus an end-to-end smoke through LiftBackStats.
// Uses the same three-exon contig as ContigTranslatorTest so lifted-coordinate math is consistent.
public class LiftBackResolverTest
{
    private static List<ContigEntry> contigMap()
    {
        return List.of(threeExonContig());
    }

    // Local unpaired builder kept on purpose: the shared fixture's primaryRecord sets paired/first-of-pair/proper-pair flags
    // which these resolver tests must not carry.
    private static SAMRecord newRecord(final String contig, final int pos, final String cigar)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        return record;
    }

    private static SAMRecord newUnmappedRecord()
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReadUnmappedFlag(true);
        return record;
    }

    @Test
    public void testUnmapped()
    {
        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(newUnmappedRecord());

        assertEquals(RecordState.UNMAPPED, result.recordState());
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testRefOnlyExonic()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "150M");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.SOLE_REF, result.decidingFeature());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(1000, result.finalPos());
        assertEquals("150M", result.finalCigar());
        assertFalse(result.hasNCigar());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxPrimaryUniqueWithinExon()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // contig pos 1, 50M -> exon 1 (chr1:100-149)

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.SOLE_TX, result.decidingFeature());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(100, result.finalPos());
        assertEquals("50M", result.finalCigar());
        assertFalse(result.hasNCigar());
    }

    @Test
    public void testTxPrimaryUniqueJunctionCrosser()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M"); // crosses exon 1 -> exon 2 (intron 200-299)

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.SOLE_TX, result.decidingFeature());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
    }

    @Test
    public void testRefTxAgree()
    {
        SAMRecord record = newRecord(CHR_1, 100, "50M"); // ref primary; Tx alt at contig pos 1 -> same locus+CIGAR
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.CONCORDANT, result.decidingFeature());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testIntronRetRefBetter()
    {
        // ref chr1:170 50M (last 20bp in intron); Tx pos 71 30M20S lifts to chr1:170-199 soft-clipped at exon boundary.
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,30M20S,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.REF_READS_THROUGH, result.decidingFeature());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxSoftClipNotAtBoundaryFallsToAmbiguous()
    {
        // ref 50M full match; Tx 25M25S - 25M ends mid-exon (chr1:194), trailing clip NOT at exon boundary -> BOTH_AMBIGUOUS.
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,25M25S,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.AMBIGUOUS, result.decidingFeature());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testMultiLocusTwoLoci()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // Tx primary + ref alt on different chrom -> two loci
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.MULTIMAPPER, result.decidingFeature());
        assertEquals(2, result.numLoci());
    }

    @Test
    public void testDistinctLocusAltBlocksMapqRescue()
    {
        // Tx primary (perfect) + a lower-scoring XA alt (3 mismatches) at a distinct locus. bwa left the read
        // MAPQ 0; two distinct genomic loci means it stays a multimapper - TARS does not override bwa's call with a
        // weaker reconstructed score - so numLoci == 2 and MAPQ is held at 0.
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setMappingQuality(0);
        record.setAttribute("XA", "chr5,+5000,50M,3;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.numLoci());
        assertEquals(0, result.updatedMapq());
    }

    @Test
    public void testNestedSpanAltCountsAsSingleLocusAndBumpsMapq()
    {
        // Regression: a 5'-softclipped isoform alt begins at a downstream exon of the SAME placement. It lifts to a
        // different start (chr1:300) but a span (300-349) nested inside the junction-crossing primary's (chr1:150-349),
        // so it is the same genomic locus, not a second one. A MAPQ-0 read here is therefore a single-locus placement
        // and must bump to 60 -- keying loci on exact start (the old behaviour) miscounted this as two loci and held it at 0.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M"); // -> chr1:150 50M100N50M, genomic span 150-349
        record.setMappingQuality(0);
        record.setAttribute("XA", TX_CONTIG + ",+101,50S50M,0;"); // -> chr1:300 50S50M, span 300-349 (nested)

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(1, result.numLoci());
        assertEquals(60, result.updatedMapq());
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

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.numLoci());
        assertEquals(0, result.updatedMapq());
    }

    @Test
    public void testRefOnlyMulti()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "50M"); // ref primary + ref alt on different chrom, no tx
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.SOLE_REF, result.decidingFeature());
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

        LiftBackResolver resolver = new LiftBackResolver(twoContigs);
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.SOLE_TX, result.decidingFeature());
        assertEquals(2, result.numLoci());
    }

    @Test
    public void testSupplementary()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "60M40H");
        record.setSupplementaryAlignmentFlag(true);

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(RecordState.SUPPLEMENTARY, result.recordState());
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.role());
    }

    @Test
    public void testSupplementaryOnTxContigGetsLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setSupplementaryAlignmentFlag(true);

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(RecordState.SUPPLEMENTARY, result.recordState());
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.role());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
    }

    @Test
    public void testSupplementaryOnTxContigUnliftablePastEnd()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M"); // pos 251 past altEnd(250)
        record.setSupplementaryAlignmentFlag(true);

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(RecordState.LIFT_FAILED, result.recordState());
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.role());
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testPrimaryUnliftablePastEnd()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(RecordState.LIFT_FAILED, result.recordState());
        assertEquals(LiftBackResult.RecordRole.PRIMARY, result.role());
        assertTrue(result.notes().contains("primary_translate_failed"));
    }

    @Test
    public void testPrimaryTrailingOverhangClampedToSoftClip()
    {
        SAMRecord record = newRecord(TX_CONTIG, 200, "100M"); // 49bp past altEnd(250) -> trailing 49S

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.SOLE_TX, result.decidingFeature());
        assertTrue(result.finalCigar().endsWith("49S"));
    }

    @Test
    public void testIntronRetRefBetterLeadingSoftClipBoundary()
    {
        // Tx alt contig pos 101 with leading 20S at exon-1/exon-2 boundary; ref chr1:300 full match.
        SAMRecord record = newRecord(CHR_1, 300, "30M");
        record.setAttribute("XA", TX_CONTIG + ",+101,20S30M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.REF_READS_THROUGH, result.decidingFeature());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testXaDedupDropsDuplicateAlts()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // two identical XA entries -> one alt retained
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,+5000,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + one deduped alt
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaDedupKeepsAltMatchingSelfButDropsXaDuplicate()
    {
        // Two XA entries lifting to the same (chr1,100,50M): one Tx, one ref. XA dedup is XA-internal only,
        // so the Tx alt is kept (drives BOTH_AGREE); the duplicate ref XA collapses.
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;" + CHR_1 + ",+100,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + Tx alt; ref XA collapsed
        assertEquals(1, result.numXaAlts());
        assertEquals(DecidingFeature.CONCORDANT, result.decidingFeature());
    }

    @Test
    public void testXaWithMalformedNmStillLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // garbled NM field must not silently drop the alt
        record.setAttribute("XA", "chr5,+5000,50M,not_a_number;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size());
        assertEquals(DecidingFeature.MULTIMAPPER, result.decidingFeature());
    }

    @Test
    public void testCrossLocusBothSplicedRemainsMultiLocus()
    {
        // Both Tx primary and ref XA alt have spliced CIGARs at different loci - ambiguous, must NOT swap.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", "chr5,+5000,50M100N50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(DecidingFeature.MULTIMAPPER, result.decidingFeature());
        LiftedAlignment self = result.liftedAlignments().stream()
                .filter(alignment -> alignment.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertTrue(self.IsPrimaryChoice);
    }

    @Test
    public void testStatsEndToEndSmoke()
    {
        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackStats stats = new LiftBackStats();

        SAMRecord r1 = newRecord(CHR_1, 1000, "150M"); // REF_SINGLE
        stats.record(r1, resolver.resolve(r1));

        SAMRecord r2 = newRecord(TX_CONTIG, 51, "100M"); // TX_SINGLE, junction-crosser
        stats.record(r2, resolver.resolve(r2));

        SAMRecord r3 = newRecord(CHR_1, 100, "50M"); // BOTH_AGREE, MAPQ=0
        r3.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");
        r3.setMappingQuality(0);
        stats.record(r3, resolver.resolve(r3));

        SAMRecord r4 = newUnmappedRecord(); // UNMAPPED
        stats.record(r4, resolver.resolve(r4));

        assertEquals(4, stats.total());
        assertEquals(1, stats.featureCount(DecidingFeature.SOLE_REF));
        assertEquals(1, stats.featureCount(DecidingFeature.SOLE_TX));
        assertEquals(1, stats.featureCount(DecidingFeature.CONCORDANT));
        assertEquals(1, stats.stateCount(RecordState.UNMAPPED));
    }

    // numLoci must reflect the deduped genomic-locus count, not the XA entry count (NH is derived from it).
    @Test
    public void testNumLociDedupesIdenticalLiftedXaEntries()
    {
        SAMRecord primary = newRecord(TX_CONTIG, 51, "100M"); // lifts to chr1:150 50M100N50M
        // four XA entries all lifting to the same locus -> numLoci still 1
        primary.setAttribute("XA",
                TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(primary);

        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testNumLociCountsDistinctLiftedLoci()
    {
        SAMRecord primary = newRecord(CHR_1, 1000, "150M");
        primary.setAttribute("XA",
                CHR_1 + ",+2000,150M,0;"
                        + CHR_1 + ",+3000,150M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(primary);

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

        LiftBackResolver noIndex = new LiftBackResolver(contigMap());
        assertEquals(0, noIndex.resolve(record).updatedMapq());

        ExonRegionIndex exonIndex = exonRegionIndex(CHR_1, List.of(new int[] { 1400, 1700 }));
        LiftBackResolver withIndex = new LiftBackResolver(contigMap(), exonIndex);
        LiftBackResult result = withIndex.resolve(record);
        assertEquals(0, result.updatedMapq());
        assertEquals(DecidingFeature.SOLE_REF, result.decidingFeature());
    }

    // Hidden tie with primary outside any annotated exon: rescue stays blocked.
    @Test
    public void testHiddenTieOutsideExonKeepsMapqZero() throws Exception
    {
        SAMRecord record = newRecord(CHR_1, 5000, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        ExonRegionIndex exonIndex = exonRegionIndex(CHR_1, List.of(new int[] { 1400, 1700 })); // exon at 1400-1700; primary at 5000 is intergenic
        LiftBackResolver resolver = new LiftBackResolver(contigMap(), exonIndex);
        assertEquals(0, resolver.resolve(record).updatedMapq());
    }

    // Direct unit tests for the extracted MAPQ policy; independent of LiftBackDiscriminator / SAMRecord plumbing.
    // decidePrimaryMapq positional args: (inputMapq, numLoci, hiddenTie, primaryFromTxContig, primaryInAnnotatedExon).
    @Test
    public void testMapqPolicy_singleLocusZeroRescues()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, false, false, false)); // single locus, MAPQ0, no hidden tie -> 60
    }

    @Test
    public void testMapqPolicy_hiddenTieRefPrimaryNoExonHoldsAtZero()
    {
        assertEquals(0, LiftBackResolver.decidePrimaryMapq(0, 1, true, false, false)); // unresolved hidden tie holds at 0
    }

    @Test
    public void testMapqPolicy_hiddenTieTxPrimaryRescues()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, true, true, false)); // tx provenance overrides hidden tie
    }

    @Test
    public void testMapqPolicy_hiddenTieInAnnotatedExonRescues()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, true, false, true)); // annotated exon overrides hidden tie
    }

    @Test
    public void testMapqPolicy_inputSixtyPassesAsRescued()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(60, 1, false, false, false)); // input 60 passes through unchanged
    }

    @Test
    public void testMapqPolicy_gradedMapqPassesThrough()
    {
        assertEquals(37, LiftBackResolver.decidePrimaryMapq(37, 1, false, false, false)); // graded signal; leave alone
    }

    @Test
    public void testMapqPolicy_multiLocusNeverBumps()
    {
        assertEquals(0, LiftBackResolver.decidePrimaryMapq(0, 2, false, false, false)); // multi-locus never bumped
    }

    @Test
    public void testMapqPolicy_randomTieNotBumped()
    {
        // a random-tie pick is a coin-flip among distinct placements, so it is not a confident unique call
        assertEquals(0, LiftBackResolver.decidePrimaryMapq(0, 1, false, false, false, true));
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, false, false, false, false));
    }

    private static LiftedAlignment liftedAt(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(
                LiftedAlignment.AlignmentSource.XA_INPUT, chrom, pos, cigar, chrom, pos, cigar,
                0, 0, null, null, null, false, true);
    }

    @Test
    public void countDistinctLociFromListRecountsPostReclaim()
    {
        // The public overload backs the emit-time NH recompute: it finds the primary via IsPrimaryChoice (as
        // buildLiftedXa does), drops Dropped alts, and collapses alts overlapping the primary - so NH stays
        // consistent with the XA after reconcileChosenPrimary's alt reclaim mutates the shared list.
        LiftedAlignment primary = liftedAt(CHR_1, 1000, "100M");   // span 1000-1099
        primary.IsPrimaryChoice = true;
        LiftedAlignment overlapping = liftedAt(CHR_1, 1050, "100M");   // 1050-1149 overlaps the primary
        LiftedAlignment distant = liftedAt(CHR_1, 5000, "100M");   // a genuinely distinct locus
        LiftedAlignment droppedDistant = liftedAt(CHR_1, 8000, "100M");
        droppedDistant.Dropped = true;

        assertEquals(1, LiftBackResolver.countDistinctLoci(List.of(primary)));
        assertEquals(1, LiftBackResolver.countDistinctLoci(List.of(primary, overlapping)));
        assertEquals(2, LiftBackResolver.countDistinctLoci(List.of(primary, distant)));
        assertEquals(1, LiftBackResolver.countDistinctLoci(List.of(primary, droppedDistant)));
        // no primary marked (e.g. an unmapped result's empty list) -> 1
        assertEquals(1, LiftBackResolver.countDistinctLoci(List.of(liftedAt(CHR_1, 1000, "100M"))));
    }

    @Test
    public void oppositeStrandXaAltsNotCollapsed()
    {
        // Two XA alts at the same genomic locus and cigar but opposite strands are distinct placements. The lifted
        // dedup key includes strand, so both survive into the XA output; previously they collapsed to one, silently
        // dropping a strand-distinct placement.
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,-5000,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        String xa = LiftBackRecordOps.buildLiftedXa(result);
        assertNotNull(xa);
        assertTrue("both strands kept: " + xa, xa.contains("chr5,+5000") && xa.contains("chr5,-5000"));
    }
}
