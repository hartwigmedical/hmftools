package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
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
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        return record;
    }

    private static SAMRecord newUnmappedRecord()
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReadUnmappedFlag(true);
        return record;
    }

    @Test
    public void testUnmapped()
    {
        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(newUnmappedRecord());

        assertEquals(LiftBackCategory.UNMAPPED, result.category());
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testRefOnlyExonic()
    {
        final SAMRecord record = newRecord(CHR_1, 1000, "150M");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.REF_SINGLE, result.category());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(1000, result.finalPos());
        assertEquals("150M", result.finalCigar());
        assertFalse(result.hasNCigar());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxPrimaryUniqueWithinExon()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // contig pos 1, 50M -> exon 1 (chr1:100-149)

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_SINGLE, result.category());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(100, result.finalPos());
        assertEquals("50M", result.finalCigar());
        assertFalse(result.hasNCigar());
    }

    @Test
    public void testTxPrimaryUniqueJunctionCrosser()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M"); // crosses exon 1 -> exon 2 (intron 200-299)

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_SINGLE, result.category());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
    }

    @Test
    public void testJunctionTxBetter()
    {
        // Tx junction-crosser primary; ref alt at same lifted locus with soft-clip.
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", CHR_1 + ",+150,50M50S,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP, result.category());
        assertTrue(result.hasNCigar());
        assertEquals(1, result.numLoci());
        // tx wins, so ref alt is dropped from the emitted set; only tx CIGAR remains at the locus.
        assertEquals(1, result.numDistinctCigarsAtPrimaryLocus());
        assertEquals(2, result.liftedAlignments().size());
        assertEquals(1, result.liftedAlignments().stream().filter(la -> !la.Dropped).count());
    }

    @Test
    public void testRefTxAgree()
    {
        final SAMRecord record = newRecord(CHR_1, 100, "50M"); // ref primary; Tx alt at contig pos 1 -> same locus+CIGAR
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_AGREE, result.category());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testIntronRetRefBetter()
    {
        // ref chr1:170 50M (last 20bp in intron); Tx pos 71 30M20S lifts to chr1:170-199 soft-clipped at exon boundary.
        final SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,30M20S,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH, result.category());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testTxSoftClipNotAtBoundaryFallsToAmbiguous()
    {
        // ref 50M full match; Tx 25M25S — 25M ends mid-exon (chr1:194), trailing clip NOT at exon boundary -> BOTH_AMBIGUOUS.
        final SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,25M25S,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_AMBIGUOUS, result.category());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testParalogMulti()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // Tx primary + ref alt on different chrom -> two loci
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI, result.category());
        assertEquals(2, result.numLoci());
    }

    @Test
    public void testSubOptimalAltDoesNotBlockMapqRescue()
    {
        // Tx primary (perfect) + sub-optimal XA alt with 3 mismatches (score 35 < 50). Sub-optimal alt is
        // not a real competitor, so numLoci == 1 and MAPQ=0 is rescued.
        final SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setMappingQuality(0);
        record.setAttribute("XA", "chr5,+5000,50M,3;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(1, result.numLoci());
        assertEquals(60, result.updatedMapq());
    }

    @Test
    public void testRefOnlyMulti()
    {
        final SAMRecord record = newRecord(CHR_1, 1000, "50M"); // ref primary + ref alt on different chrom, no tx
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.REF_MULTI, result.category());
        assertEquals(2, result.numLoci());
    }

    @Test
    public void testTxPrimaryMulti()
    {
        final ContigEntry entryB = new ContigEntry(
                "ensG_OTHER_T", 1, 100, "GO", "OTHER", "TO", "chr5", 1,
                List.of(new BaseRegion(2000, 2099)));
        final List<ContigEntry> twoContigs = List.of(threeExonContig(), entryB);

        final SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "ensG_OTHER_T,+1,50M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(twoContigs);
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_MULTI, result.category());
        assertEquals(2, result.numLoci());
    }

    @Test
    public void testSupplementary()
    {
        final SAMRecord record = newRecord(CHR_1, 1000, "60M40H");
        record.setSupplementaryAlignmentFlag(true);

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.SUPPLEMENTARY, result.category());
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.role());
    }

    @Test
    public void testSupplementaryOnTxContigGetsLifted()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setSupplementaryAlignmentFlag(true);

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.SUPPLEMENTARY, result.category());
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.role());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
    }

    @Test
    public void testSupplementaryOnTxContigUnliftablePastEnd()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 251, "10M"); // pos 251 past altEnd(250)
        record.setSupplementaryAlignmentFlag(true);

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.LIFT_FAILED, result.category());
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.role());
        assertTrue(result.liftedAlignments().isEmpty());
    }

    @Test
    public void testPrimaryUnliftablePastEnd()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 251, "10M");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.LIFT_FAILED, result.category());
        assertEquals(LiftBackResult.RecordRole.PRIMARY, result.role());
        assertTrue(result.notes().contains("primary_translate_failed"));
    }

    @Test
    public void testPrimaryTrailingOverhangClampedToSoftClip()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 200, "100M"); // 49bp past altEnd(250) -> trailing 49S

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_SINGLE, result.category());
        assertTrue(result.finalCigar().endsWith("49S"));
    }

    @Test
    public void testIntronRetRefBetterLeadingSoftClipBoundary()
    {
        // Tx alt contig pos 101 with leading 20S at exon-1/exon-2 boundary; ref chr1:300 full match.
        final SAMRecord record = newRecord(CHR_1, 300, "30M");
        record.setAttribute("XA", TX_CONTIG + ",+101,20S30M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH, result.category());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testXaDedupDropsDuplicateAlts()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // two identical XA entries -> one alt retained
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,+5000,50M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + one deduped alt
        assertEquals(1, result.numXaAlts());
    }

    @Test
    public void testXaDedupKeepsAltMatchingSelfButDropsXaDuplicate()
    {
        // Two XA entries lifting to the same (chr1,100,50M): one Tx, one ref. XA dedup is XA-internal only,
        // so the Tx alt is kept (drives BOTH_AGREE); the duplicate ref XA collapses.
        final SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;" + CHR_1 + ",+100,50M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size()); // self + Tx alt; ref XA collapsed
        assertEquals(1, result.numXaAlts());
        assertEquals(LiftBackCategory.BOTH_AGREE, result.category());
    }

    @Test
    public void testXaWithMalformedNmStillLifted()
    {
        final SAMRecord record = newRecord(TX_CONTIG, 1, "50M"); // garbled NM field must not silently drop the alt
        record.setAttribute("XA", "chr5,+5000,50M,not_a_number;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.liftedAlignments().size());
        assertEquals(LiftBackCategory.BOTH_MULTI, result.category());
    }

    @Test
    public void testCrossLocusFavoursTxSwapsToTx()
    {
        // bwa picked paralog chr5:5000 as primary; Tx XA lifts to chr1:150 50M100N50M. Expect swap to Tx, MAPQ rescued.
        final SAMRecord record = newRecord("chr5", 5000, "100M");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;");
        record.setMappingQuality(60);

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI_TX_JUNCTION, result.category());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
        assertEquals(60, result.updatedMapq());
        assertTrue(result.notes().contains("swapped_ref_to_tx"));

        final LiftedAlignment originalSelf = result.liftedAlignments().stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertFalse(originalSelf.IsPrimaryChoice);
        assertFalse(originalSelf.Dropped);

        final LiftedAlignment winner = result.liftedAlignments().stream()
                .filter(la -> la.IsPrimaryChoice)
                .findFirst().orElseThrow();
        assertTrue(winner.fromTxContig());
        assertEquals("50M100N50M", winner.LiftedCigar);
    }

    @Test
    public void testCrossLocusBothSplicedRemainsMultiLocus()
    {
        // Both Tx primary and ref XA alt have spliced CIGARs at different loci — ambiguous, must NOT swap.
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", "chr5,+5000,50M100N50M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI, result.category());
        final LiftedAlignment self = result.liftedAlignments().stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertTrue(self.IsPrimaryChoice);
    }

    @Test
    public void testCrossLocusFavoursTxDropsOtherRefAlts()
    {
        // Tx XA wins (chr1:150 spliced); original self (chr5) kept as informative paralog; chr7 ref alt dropped.
        final SAMRecord record = newRecord("chr5", 5000, "100M");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;chr7,+7000,100M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI_TX_JUNCTION, result.category());

        final LiftedAlignment chr7Alt = result.liftedAlignments().stream()
                .filter(la -> "chr7".equals(la.LiftedChrom))
                .findFirst().orElseThrow();
        assertTrue(chr7Alt.Dropped);

        final LiftedAlignment chr5Self = result.liftedAlignments().stream()
                .filter(la -> "chr5".equals(la.LiftedChrom))
                .findFirst().orElseThrow();
        assertFalse(chr5Self.Dropped);
        assertFalse(chr5Self.IsPrimaryChoice);
    }

    @Test
    public void testJunctionFavoursTxSwapsWhenSelfIsRef()
    {
        // bwa picked ref soft-clipped (50M50S) as primary; Tx XA lifts to the same locus with 50M100N50M. Tx wins and swaps.
        final SAMRecord record = newRecord(CHR_1, 150, "50M50S");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP, result.category());
        assertEquals("50M100N50M", result.finalCigar());
        assertTrue(result.hasNCigar());
        assertTrue(result.notes().contains("swapped_ref_to_tx"));
    }

    @Test
    public void testJunctionRefMatchKeepsRefAndDropsTx()
    {
        // Ref full-match 151M through the supposed intron is overwhelming evidence of an unspliced read; keep ref, drop Tx.
        final SAMRecord record = newRecord(CHR_1, 150, "151M");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH, result.category());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("151M", result.finalCigar());
        assertFalse(result.hasNCigar());

        final LiftedAlignment self = result.liftedAlignments().stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertTrue(self.IsPrimaryChoice);
        assertFalse(self.Dropped);

        final LiftedAlignment txAlt = result.liftedAlignments().stream()
                .filter(LiftedAlignment::fromTxContig)
                .findFirst().orElseThrow();
        assertTrue(txAlt.Dropped);
        assertFalse(txAlt.IsPrimaryChoice);
    }

    @Test
    public void testJunctionRefMatchSwapsToRefWhenSelfIsTx()
    {
        // Tx primary; ref XA is unspliced full-match. Discriminator swaps to ref and drops Tx alt.
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", CHR_1 + ",+150,151M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH, result.category());
        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("151M", result.finalCigar());
        assertFalse(result.hasNCigar());
        assertTrue(result.notes().contains("swapped_tx_to_ref"));

        final LiftedAlignment txSelf = result.liftedAlignments().stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertFalse(txSelf.IsPrimaryChoice);
        assertFalse(txSelf.Dropped); // demoted to informative XA, not dropped

        final LiftedAlignment winner = result.liftedAlignments().stream()
                .filter(la -> la.IsPrimaryChoice)
                .findFirst().orElseThrow();
        assertFalse(winner.fromTxContig());
    }

    @Test
    public void testStatsEndToEndSmoke()
    {
        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackStats stats = new LiftBackStats();

        final SAMRecord r1 = newRecord(CHR_1, 1000, "150M"); // REF_SINGLE
        stats.record(r1, resolver.resolve(r1));

        final SAMRecord r2 = newRecord(TX_CONTIG, 51, "100M"); // TX_SINGLE, junction-crosser
        stats.record(r2, resolver.resolve(r2));

        final SAMRecord r3 = newRecord(CHR_1, 100, "50M"); // BOTH_AGREE, MAPQ=0
        r3.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");
        r3.setMappingQuality(0);
        stats.record(r3, resolver.resolve(r3));

        final SAMRecord r4 = newUnmappedRecord(); // UNMAPPED
        stats.record(r4, resolver.resolve(r4));

        assertEquals(4, stats.total());
        assertEquals(1, stats.categoryCount(LiftBackCategory.REF_SINGLE));
        assertEquals(1, stats.categoryCount(LiftBackCategory.TX_SINGLE));
        assertEquals(1, stats.categoryCount(LiftBackCategory.BOTH_AGREE));
        assertEquals(1, stats.categoryCount(LiftBackCategory.UNMAPPED));
    }

    // numLoci must reflect the deduped genomic-locus count, not the XA entry count (NH is derived from it).
    @Test
    public void testNumLociDedupesIdenticalLiftedXaEntries()
    {
        final SAMRecord primary = newRecord(TX_CONTIG, 51, "100M"); // lifts to chr1:150 50M100N50M
        // four XA entries all lifting to the same locus -> numLoci still 1
        primary.setAttribute("XA",
                TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;"
                        + TX_CONTIG + ",+51,100M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(primary);

        assertEquals(CHR_1, result.finalChrom());
        assertEquals(150, result.finalPos());
        assertEquals("50M100N50M", result.finalCigar());
        assertEquals(1, result.numLoci());
    }

    @Test
    public void testNumLociCountsDistinctLiftedLoci()
    {
        final SAMRecord primary = newRecord(CHR_1, 1000, "150M");
        primary.setAttribute("XA",
                CHR_1 + ",+2000,150M,0;"
                        + CHR_1 + ",+3000,150M,0;");

        final LiftBackResolver resolver = new LiftBackResolver(contigMap());
        final LiftBackResult result = resolver.resolve(primary);

        assertEquals(3, result.numLoci());
    }

    // Builds an ExonRegionIndex from an in-memory ensembl cache (shared by rescue-path tests).
    private static ExonRegionIndex buildExonIndex(final List<int[]> exons)
    {
        final EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, CHR_1, List.of(createEnsemblGeneData("ENSG_TEST", "TESTG", CHR_1, 1, 1, 100000)));

        final TranscriptData transcript = new TranscriptData(
                1, "ENST_TEST", "ENSG_TEST", true, (byte) 1, 1, 100000, null, null, "protein_coding", "");
        final List<ExonData> exonData = new ArrayList<>();
        int rank = 1;
        for(final int[] exon : exons)
            exonData.add(new ExonData(1, exon[0], exon[1], rank++, -1, -1));
        transcript.setExons(exonData);
        addTransExonData(cache, "ENSG_TEST", List.of(transcript));

        return ExonRegionIndex.fromCache(cache, V38);
    }

    // Hidden tie (XS==AS) on a ref-only primary: rescue is gated on a tx match, so exon evidence cannot bump MAPQ.
    @Test
    public void testHiddenTieRefOnlyNeverRescuesMapq() throws Exception
    {
        final SAMRecord record = newRecord(CHR_1, 1500, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        final LiftBackResolver noIndex = new LiftBackResolver(contigMap());
        assertEquals(0, noIndex.resolve(record).updatedMapq());

        final ExonRegionIndex idx = buildExonIndex(List.of(new int[] { 1400, 1700 })); // covers pos 1500
        final LiftBackResolver withIndex = new LiftBackResolver(contigMap(), idx);
        final LiftBackResult result = withIndex.resolve(record);
        assertEquals(0, result.updatedMapq());
        assertEquals(LiftBackCategory.REF_SINGLE, result.category());
    }

    // Hidden tie with primary outside any annotated exon: rescue stays blocked.
    @Test
    public void testHiddenTieOutsideExonKeepsMapqZero() throws Exception
    {
        final SAMRecord record = newRecord(CHR_1, 5000, "150M");
        record.setMappingQuality(0);
        record.setAttribute("AS", 151);
        record.setAttribute("XS", 151);

        final ExonRegionIndex idx = buildExonIndex(List.of(new int[] { 1400, 1700 })); // exon at 1400-1700; primary at 5000 is intergenic
        final LiftBackResolver resolver = new LiftBackResolver(contigMap(), idx);
        assertEquals(0, resolver.resolve(record).updatedMapq());
    }

    // Direct unit tests for the extracted MAPQ policy; independent of LiftBackDiscriminator / SAMRecord plumbing.
    // decidePrimaryMapq positional args: (inputMapq, numLoci, swapped, hiddenTie, primaryFromTxContig, primaryInAnnotatedExon, hasTxMatch).
    @Test
    public void testMapqPolicy_swappedAlwaysRescues()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, true, true, false, false, true)); // swap is decisive even with hidden tie
    }

    @Test
    public void testMapqPolicy_singleLocusZeroRescues()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, false, false, false, false, true));
    }

    @Test
    public void testMapqPolicy_hiddenTieRefPrimaryNoExonHoldsAtZero()
    {
        assertEquals(0, LiftBackResolver.decidePrimaryMapq(0, 1, false, true, false, false, true)); // unresolved hidden tie
    }

    @Test
    public void testMapqPolicy_hiddenTieTxPrimaryRescues()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, false, true, true, false, true)); // tx provenance overrides hidden tie
    }

    @Test
    public void testMapqPolicy_hiddenTieInAnnotatedExonRescues()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(0, 1, false, true, false, true, true)); // exon evidence overrides hidden tie
    }

    @Test
    public void testMapqPolicy_inputSixtyPassesAsRescued()
    {
        assertEquals(60, LiftBackResolver.decidePrimaryMapq(60, 1, false, false, false, false, true)); // no-op since RESCUED_MAPQ==60
    }

    @Test
    public void testMapqPolicy_gradedMapqPassesThrough()
    {
        assertEquals(37, LiftBackResolver.decidePrimaryMapq(37, 1, false, false, false, false, true)); // graded signal; leave alone
    }

    @Test
    public void testMapqPolicy_multiLocusZeroStaysZero()
    {
        assertEquals(0, LiftBackResolver.decidePrimaryMapq(0, 2, false, false, false, false, true)); // real alt exists; honest 0 stands
    }

    @Test
    public void testMapqPolicy_noTxMatchNeverBumps()
    {
        assertEquals(0, LiftBackResolver.decidePrimaryMapq(0, 1, false, false, false, false, false)); // no tx -> MAPQ-0 is honest
        assertEquals(0, LiftBackResolver.decidePrimaryMapq(0, 1, true, false, false, false, false)); // swap can't rescue without tx match
    }
}
