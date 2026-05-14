package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

// per-category sanity checks on LiftBackResolver + a small end-to-end smoke through LiftBackStats.
// uses the same three-exon test contig as ContigTranslatorTest so lifted-coordinate math is consistent.
public class LiftBackResolverTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    // exon spans on chr1: 100-199, 300-399, 500-549. introns: 200-299, 400-499. contig length 250.
    // alt contig has the transcript laid down at positions 1..250 so test positions read like local coords.
    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    private static List<ContigEntry> contigMap()
    {
        return List.of(threeExonContig());
    }

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

        assertEquals(LiftBackCategory.UNMAPPED, result.Category);
        assertTrue(result.LiftedAlignments.isEmpty());
    }

    @Test
    public void testRefOnlyExonic()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "150M");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.REF_SINGLE, result.Category);
        assertEquals(CHR_1, result.FinalChrom);
        assertEquals(1000, result.FinalPos);
        assertEquals("150M", result.FinalCigar);
        assertFalse(result.HasNCigar);
        assertEquals(1, result.NumLoci);
    }

    @Test
    public void testTxPrimaryUniqueWithinExon()
    {
        // contig pos 1, 50M — entirely within exon 1 (chr1:100-149)
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_SINGLE, result.Category);
        assertEquals(CHR_1, result.FinalChrom);
        assertEquals(100, result.FinalPos);
        assertEquals("50M", result.FinalCigar);
        assertFalse(result.HasNCigar);
    }

    @Test
    public void testTxPrimaryUniqueJunctionCrosser()
    {
        // contig pos 51, 100M — crosses exon 1 -> exon 2 (intron 200-299, length 100)
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_SINGLE, result.Category);
        assertEquals(CHR_1, result.FinalChrom);
        assertEquals(150, result.FinalPos);
        assertEquals("50M100N50M", result.FinalCigar);
        assertTrue(result.HasNCigar);
    }

    @Test
    public void testJunctionTxBetter()
    {
        // primary on Tx (junction-crosser), ref alt in XA at same lifted locus with soft-clip
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", CHR_1 + ",+150,50M50S,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP, result.Category);
        assertTrue(result.HasNCigar);
        assertEquals(1, result.NumLoci);
        // ref alt is dropped from the emitted set because tx wins this discriminator outcome,
        // so only the tx CIGAR remains at the primary locus. The dropped ref alt stays in
        // LiftedAlignments for TSV-B diagnostics.
        assertEquals(1, result.NumDistinctCigarsAtPrimaryLocus);
        assertEquals(2, result.LiftedAlignments.size());
        assertEquals(1, result.LiftedAlignments.stream().filter(la -> !la.Dropped).count());
    }

    @Test
    public void testRefTxAgree()
    {
        // primary on ref at chr1:100, Tx alt at contig pos 1 -> same lifted locus, same CIGAR
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_AGREE, result.Category);
        assertEquals(1, result.NumLoci);
    }

    @Test
    public void testIntronRetRefBetter()
    {
        // primary on ref full match overlapping exon 1 + intron 1; Tx alt soft-clipped at exon boundary
        // ref: chr1:170, 50M -> covers chr1:170-219 (last 20 in intron)
        // Tx: contig pos 71, 30M20S -> 30M lifts to chr1:170-199 (end of exon 1), 20S past contig boundary
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,30M20S,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH, result.Category);
        assertEquals(1, result.NumLoci);
    }

    @Test
    public void testTxSoftClipNotAtBoundaryFallsToAmbiguous()
    {
        // ref pos 170, 50M (chr1:170-219, full match). Tx pos 71, 25M25S — 25M lifts to chr1:170-194, ending
        // mid exon-1 (well shy of exon-1 end 199), so trailing 25S is NOT at any interior exon boundary.
        // Same lifted locus (chr1:170), distinct CIGARs, refFullMatch=true, txSoftClipAtBoundary=false -> BOTH_AMBIGUOUS.
        SAMRecord record = newRecord(CHR_1, 170, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+71,25M25S,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_AMBIGUOUS, result.Category);
        assertEquals(1, result.NumLoci);
    }

    @Test
    public void testParalogMulti()
    {
        // primary on Tx, ref alt at a different chromosome -> two distinct loci, ref present
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI, result.Category);
        assertEquals(2, result.NumLoci);
    }

    @Test
    public void testRefOnlyMulti()
    {
        // primary on ref, XA alt on a different ref chromosome -> two distinct loci, no tx
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.REF_MULTI, result.Category);
        assertEquals(2, result.NumLoci);
    }

    @Test
    public void testTxPrimaryMulti()
    {
        // primary on Tx, alt on a different ref locus, no actual ref alt -> Tx-only multi-locus
        // (using two distinct Tx contigs would need a second ContigEntry; instead use a different contig name
        //  not in our map so it's treated as ref — defeats the test. Simulate with two different lifted positions
        //  of two Tx alts. Since our map only has one contig, instead build a multi-contig test:)
        ContigEntry entryB = new ContigEntry(
                "ensG_OTHER_T", 1, 100, "GO", "OTHER", "TO", "chr5",
                List.of(new BaseRegion(2000, 2099)));
        List<ContigEntry> twoContigs = List.of(threeExonContig(), entryB);

        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "ensG_OTHER_T,+1,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(twoContigs);
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_MULTI, result.Category);
        assertEquals(2, result.NumLoci);
    }

    @Test
    public void testSupplementary()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "60M40H");
        record.setSupplementaryAlignmentFlag(true);

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.SUPPLEMENTARY, result.Category);
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.Role);
    }

    @Test
    public void testSupplementaryOnTxContigGetsLifted()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setSupplementaryAlignmentFlag(true);

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.SUPPLEMENTARY, result.Category);
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.Role);
        assertEquals(CHR_1, result.FinalChrom);
        assertEquals(150, result.FinalPos);
        assertEquals("50M100N50M", result.FinalCigar);
        assertTrue(result.HasNCigar);
    }

    @Test
    public void testSupplementaryOnTxContigUnliftablePastEnd()
    {
        // pos 251 is past altEnd (250) entirely — no overhang clamp can save it
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");
        record.setSupplementaryAlignmentFlag(true);

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.LIFT_FAILED, result.Category);
        assertEquals(LiftBackResult.RecordRole.SUPPLEMENTARY, result.Role);
        assertTrue(result.LiftedAlignments.isEmpty());
    }

    @Test
    public void testPrimaryUnliftablePastEnd()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.LIFT_FAILED, result.Category);
        assertEquals(LiftBackResult.RecordRole.PRIMARY, result.Role);
        assertTrue(result.Notes.contains("primary_translate_failed"));
    }

    @Test
    public void testPrimaryTrailingOverhangClampedToSoftClip()
    {
        // pos 200 + 100M extends 49 bases past altEnd (250) — clamp converts overhang to trailing S
        SAMRecord record = newRecord(TX_CONTIG, 200, "100M");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.TX_SINGLE, result.Category);
        assertTrue(result.FinalCigar.endsWith("49S"));
    }

    @Test
    public void testIntronRetRefBetterLeadingSoftClipBoundary()
    {
        // Tx alt at exon-2 start (contig pos 101) with leading 20S abutting the exon-1/exon-2 boundary.
        // ref alt at chr1:300 full match -> single locus, refFullMatch + txSoftClipAtBoundary -> BOTH_TX_SOFTCLIP_REF_MATCH.
        SAMRecord record = newRecord(CHR_1, 300, "30M");
        record.setAttribute("XA", TX_CONTIG + ",+101,20S30M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH, result.Category);
        assertEquals(1, result.NumLoci);
    }

    @Test
    public void testXaDedupDropsDuplicateAlts()
    {
        // two identical XA entries -> only one lifted alt retained
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,0;chr5,+5000,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        // self + one deduped alt = 2 alignments
        assertEquals(2, result.LiftedAlignments.size());
        assertEquals(1, result.NumXaAlts);
    }

    @Test
    public void testXaDedupKeepsAltMatchingSelfButDropsXaDuplicate()
    {
        // self at chr1:100/50M; two XA entries that both lift to the same (chr1, 100, 50M) — one Tx, one ref.
        // Spec: dedup XA among itself only, NOT against self. The Tx alt is kept (so BOTH_AGREE can fire),
        // the duplicate ref XA is collapsed.
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("XA", TX_CONTIG + ",+1,50M,0;" + CHR_1 + ",+100,50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        // self + 1 surviving XA alt (the Tx one; the chr1,100,50M XA collapses to the Tx one as both share lifted key)
        assertEquals(2, result.LiftedAlignments.size());
        assertEquals(1, result.NumXaAlts);
        assertEquals(LiftBackCategory.BOTH_AGREE, result.Category);
    }

    @Test
    public void testXaWithMalformedNmStillLifted()
    {
        // bwa XA NM field is sometimes garbled; alt should still be lifted, not silently dropped
        SAMRecord record = newRecord(TX_CONTIG, 1, "50M");
        record.setAttribute("XA", "chr5,+5000,50M,not_a_number;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(2, result.LiftedAlignments.size());
        assertEquals(LiftBackCategory.BOTH_MULTI, result.Category);
    }

    @Test
    public void testCrossLocusFavoursTxSwapsToTx()
    {
        // bwa picked an intronless paralog (ref locus chr5:5000, full match) as primary; tx XA lifts to
        // chr1:150 with the spliced CIGAR 50M100N50M. Expect: category CROSS_LOCUS_FAVOURS_TX, BAM primary
        // swapped to the tx locus (chr1:150 with N), original ref demoted to XA, MAPQ rescued.
        SAMRecord record = newRecord("chr5", 5000, "100M");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;");
        record.setMappingQuality(60);

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI_TX_JUNCTION, result.Category);
        assertEquals(CHR_1, result.FinalChrom);
        assertEquals(150, result.FinalPos);
        assertEquals("50M100N50M", result.FinalCigar);
        assertTrue(result.HasNCigar);
        assertEquals(60, result.UpdatedMapq);
        assertTrue(result.Notes.contains("swapped_ref_to_tx"));

        // self (chr5 ref) should still be in the lifted set but not the primary choice and not dropped
        LiftedAlignment originalSelf = result.LiftedAlignments.stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertFalse(originalSelf.IsPrimaryChoice);
        assertFalse(originalSelf.Dropped);

        // the tx XA alt is now the primary choice
        LiftedAlignment winner = result.LiftedAlignments.stream()
                .filter(la -> la.IsPrimaryChoice)
                .findFirst().orElseThrow();
        assertTrue(winner.fromTxContig());
        assertEquals("50M100N50M", winner.LiftedCigar);
    }

    @Test
    public void testCrossLocusBothSplicedRemainsMultiLocus()
    {
        // primary on tx (junction-crosser), XA alt is a ref alignment that also has an N-CIGAR at a
        // different locus. Both sides found splices — ambiguous, must NOT swap. Falls into MULTI_LOCUS.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", "chr5,+5000,50M100N50M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI, result.Category);
        // self stays primary (no swap)
        LiftedAlignment self = result.LiftedAlignments.stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertTrue(self.IsPrimaryChoice);
    }

    @Test
    public void testCrossLocusFavoursTxDropsOtherRefAlts()
    {
        // self ref at chr5 (paralog A), tx XA at chr1:150 spliced (winner), another ref XA at chr7
        // (paralog B). Expect winner promoted; original self goes to XA (informative paralog); chr7 ref
        // alt gets Dropped (paralog noise).
        SAMRecord record = newRecord("chr5", 5000, "100M");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;chr7,+7000,100M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_MULTI_TX_JUNCTION, result.Category);

        LiftedAlignment chr7Alt = result.LiftedAlignments.stream()
                .filter(la -> "chr7".equals(la.LiftedChrom))
                .findFirst().orElseThrow();
        assertTrue(chr7Alt.Dropped);

        LiftedAlignment chr5Self = result.LiftedAlignments.stream()
                .filter(la -> "chr5".equals(la.LiftedChrom))
                .findFirst().orElseThrow();
        assertFalse(chr5Self.Dropped);
        assertFalse(chr5Self.IsPrimaryChoice);
    }

    @Test
    public void testJunctionFavoursTxSwapsWhenSelfIsRef()
    {
        // single-locus junction case but bwa picked ref as primary: ref self soft-clipped at chr1:150 with
        // 50M50S, tx XA at contig pos 51 lifts to chr1:150 with 50M100N50M. Tx wins the discriminator and
        // the BAM primary now actually swaps (previously emitted "self_was_ref_no_swap" without swapping).
        SAMRecord record = newRecord(CHR_1, 150, "50M50S");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP, result.Category);
        assertEquals("50M100N50M", result.FinalCigar);
        assertTrue(result.HasNCigar);
        assertTrue(result.Notes.contains("swapped_ref_to_tx"));
    }

    @Test
    public void testJunctionRefMatchKeepsRefAndDropsTx()
    {
        // single-locus: ref self full-match 151M at chr1:150, tx XA at contig pos 51 lifts to 50M100N50M.
        // Ref full-match through the supposed intron is overwhelming evidence the read is unspliced —
        // we keep ref and drop the tx alt (opposite of BOTH_TX_JUNCTION_REF_SOFTCLIP).
        SAMRecord record = newRecord(CHR_1, 150, "151M");
        record.setAttribute("XA", TX_CONTIG + ",+51,100M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH, result.Category);
        assertEquals(CHR_1, result.FinalChrom);
        assertEquals(150, result.FinalPos);
        assertEquals("151M", result.FinalCigar);
        assertFalse(result.HasNCigar);

        // ref self stays primary
        LiftedAlignment self = result.LiftedAlignments.stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertTrue(self.IsPrimaryChoice);
        assertFalse(self.Dropped);

        // tx alt is dropped
        LiftedAlignment txAlt = result.LiftedAlignments.stream()
                .filter(LiftedAlignment::fromTxContig)
                .findFirst().orElseThrow();
        assertTrue(txAlt.Dropped);
        assertFalse(txAlt.IsPrimaryChoice);
    }

    @Test
    public void testJunctionRefMatchSwapsToRefWhenSelfIsTx()
    {
        // bwa primary is the tx alignment (51M+99M lifted to 50M100N50M), ref XA is the unspliced full-match.
        // The discriminator should swap onto the ref locus and drop the tx alt.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("XA", CHR_1 + ",+150,151M,0;");

        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackResult result = resolver.resolve(record);

        assertEquals(LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH, result.Category);
        assertEquals(CHR_1, result.FinalChrom);
        assertEquals(150, result.FinalPos);
        assertEquals("151M", result.FinalCigar);
        assertFalse(result.HasNCigar);
        assertTrue(result.Notes.contains("swapped_tx_to_ref"));

        // tx self gets demoted (kept in XA, not dropped — informative tx hit)
        LiftedAlignment txSelf = result.LiftedAlignments.stream()
                .filter(la -> la.Source == LiftedAlignment.AlignmentSource.SELF)
                .findFirst().orElseThrow();
        assertFalse(txSelf.IsPrimaryChoice);
        assertFalse(txSelf.Dropped);

        // ref alt is now the primary
        LiftedAlignment winner = result.LiftedAlignments.stream()
                .filter(la -> la.IsPrimaryChoice)
                .findFirst().orElseThrow();
        assertFalse(winner.fromTxContig());
    }

    @Test
    public void testStatsEndToEndSmoke()
    {
        LiftBackResolver resolver = new LiftBackResolver(contigMap());
        LiftBackStats stats = new LiftBackStats();

        // r1: ref-only exonic, MAPQ=60, no XA -> REF_SINGLE, REF_ONLY, MAPQ_POS_UNIQUE
        SAMRecord r1 = newRecord(CHR_1, 1000, "150M");
        stats.record(r1, resolver.resolve(r1));

        // r2: Tx primary junction-crosser, MAPQ=60, no XA -> TX_SINGLE, TX_ONLY, MAPQ_POS_UNIQUE
        SAMRecord r2 = newRecord(TX_CONTIG, 51, "100M");
        stats.record(r2, resolver.resolve(r2));

        // r3: ref + Tx agree, MAPQ=0 -> BOTH_AGREE, REF_AND_TX, MAPQ_ZERO
        SAMRecord r3 = newRecord(CHR_1, 100, "50M");
        r3.setAttribute("XA", TX_CONTIG + ",+1,50M,0;");
        r3.setMappingQuality(0);
        stats.record(r3, resolver.resolve(r3));

        // r4: unmapped -> UNMAPPED, NONE composition
        SAMRecord r4 = newUnmappedRecord();
        stats.record(r4, resolver.resolve(r4));

        assertEquals(4, stats.total());
        assertEquals(1, stats.categoryCount(LiftBackCategory.REF_SINGLE));
        assertEquals(1, stats.categoryCount(LiftBackCategory.TX_SINGLE));
        assertEquals(1, stats.categoryCount(LiftBackCategory.BOTH_AGREE));
        assertEquals(1, stats.categoryCount(LiftBackCategory.UNMAPPED));
    }

}
