package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.splice.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBack.applyResultToRecord;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBack.toLiftedMateInfo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

// black-box end-to-end harness around the full SpliceLiftBack record pipeline (resolver + apply + SA
// rewrite + mate patch). A test prepares a fragment (R1, optional R2, optional supplementary R1), hands
// it to processFragment(), then asserts on the final emitted SAMRecord fields. No inner-class poking,
// no resolver wiring — write scenarios in terms of "what bwa produced" and "what should land in the
// output BAM".
public class SpliceLiftBackEndToEndTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;
    private static final String XA_TAG = "XA";

    // exon spans on chr1: 100-199, 300-399, 500-549. introns: 200-299, 400-499. contig len 250.
    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    private static SAMRecord newPrimary(
            final String readName, final boolean firstOfPair, final String contig, final int pos, final String cigar)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(firstOfPair);
        record.setProperPairFlag(true);
        return record;
    }

    private static SAMRecord newSupplementary(
            final String readName, final boolean firstOfPair, final String contig, final int pos, final String cigar)
    {
        final SAMRecord record = newPrimary(readName, firstOfPair, contig, pos, cigar);
        record.setSupplementaryAlignmentFlag(true);
        return record;
    }

    // runs the full per-record pipeline (resolver + apply + SA rewrite + mate patch) over a list of
    // records belonging to the SAME read pair, mirroring what pass 1 + pass 2 of SpliceLiftBack do.
    // Mutates the input records in place and returns them in input order.
    private static List<SAMRecord> processFragment(final List<ContigEntry> entries, final List<SAMRecord> records)
    {
        final LiftBackResolver resolver = new LiftBackResolver(entries);
        final LiftedMateInfoCache cache = new LiftedMateInfoCache();

        // pass 1: cache lifted primary info per paired primary
        for(final SAMRecord record : records)
        {
            if(!record.getReadPairedFlag() || record.isSecondaryOrSupplementary())
                continue;
            final LiftBackResult result = resolver.resolve(record);
            cache.recordPrimaryAlignment(record.getReadName(), record.getFirstOfPairFlag(), toLiftedMateInfo(record, result));
        }

        // pass 2: lift, rewrite SA, patch mate fields
        for(final SAMRecord record : records)
        {
            final LiftBackResult result = resolver.resolve(record);
            applyResultToRecord(record, result, cache);
            record.setAttribute(SA_ATTRIBUTE, rewriteSaTag(record.getStringAttribute(SA_ATTRIBUTE), resolver));
            patchMateFields(record, cache);
        }

        return records;
    }

    @Test
    public void txPrimaryReadsBecomeGenomicWithJunctionCigar()
    {
        // single R1, primary on tx contig spanning exon1 -> exon2
        final SAMRecord r1 = newPrimary("read1", true, TX_CONTIG, 51, "100M");
        // R2 on the same tx contig spanning exon2 tail -> exon3. Start at tx 197 so the leading
        // exon2 anchor is 4M (above ContigTranslator.trimMicroAnchors threshold of 3); tx 200
        // would lift to a 1M head that the trim would drop.
        final SAMRecord r2 = newPrimary("read1", false, TX_CONTIG, 197, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        // R1: junction-crossing CIGAR
        assertEquals(CHR_1, r1.getReferenceName());
        assertEquals(150, r1.getAlignmentStart());
        assertEquals("50M100N50M", r1.getCigarString());

        // R1 mate fields refer to R2's lifted coords on chr1
        assertEquals(CHR_1, r1.getMateReferenceName());
        assertEquals(396, r1.getMateAlignmentStart());

        // R2 lifted into exon2/intron2/exon3 region with its own junction
        assertEquals(CHR_1, r2.getReferenceName());
        assertEquals(396, r2.getAlignmentStart());
        assertTrue(r2.getCigarString().contains("N"));

        // R2 mate fields refer back to R1
        assertEquals(150, r2.getMateAlignmentStart());
    }

    @Test
    public void supplementarySaTagRewrittenToGenomicCoords()
    {
        // R1 primary on tx contig with a supplementary on the same tx contig further along. The
        // supp starts at tx 197 so it has a 4M anchor in exon2 before the junction — above
        // ContigTranslator.trimMicroAnchors threshold of 3.
        final SAMRecord r1 = newPrimary("read2", true, TX_CONTIG, 51, "100M");
        r1.setAttribute(SA_ATTRIBUTE, TX_CONTIG + ",197,+,50M,60,0;");

        final SAMRecord r1Supp = newSupplementary("read2", true, TX_CONTIG, 197, "50M");
        r1Supp.setAttribute(SA_ATTRIBUTE, TX_CONTIG + ",51,+,100M,60,0;");

        final SAMRecord r2 = newPrimary("read2", false, CHR_1, 600, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r1Supp, r2));

        // r1's SA entry pointed to (TX_CONTIG, 197, +, 50M) which lifts to (chr1, 396, +, ...) per the exon spans
        final String r1Sa = r1.getStringAttribute(SA_ATTRIBUTE);
        assertTrue("expected lifted SA entry on chr1, got: " + r1Sa, r1Sa != null && r1Sa.startsWith(CHR_1 + ","));
        // no _tx contig name should leak through
        assertFalse(r1Sa.contains(TX_CONTIG));

        // supplementary primary itself is lifted to chr1 too
        assertEquals(CHR_1, r1Supp.getReferenceName());
        assertEquals(396, r1Supp.getAlignmentStart());

        // r1 and r1Supp share the same mate (r2 on chr1:600), so both carry consistent mate info
        assertEquals(CHR_1, r1.getMateReferenceName());
        assertEquals(600, r1.getMateAlignmentStart());
        assertEquals(CHR_1, r1Supp.getMateReferenceName());
        assertEquals(600, r1Supp.getMateAlignmentStart());
    }

    @Test
    public void unliftablePrimaryFlagsRecordUnmappedAndClearsMatesProperPair()
    {
        // R1 starts past altEnd — translation fails -> emit as unmapped
        final SAMRecord r1 = newPrimary("read3", true, TX_CONTIG, 251, "50M");
        // R2 lifts cleanly on the genome
        final SAMRecord r2 = newPrimary("read3", false, CHR_1, 600, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertTrue(r1.getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, r1.getCigarString());
        assertNull(r1.getStringAttribute(XA_TAG));

        // R2 stays mapped, but its mate (R1) is unmapped -> proper-pair must be cleared
        assertEquals(CHR_1, r2.getReferenceName());
        assertFalse(r2.getProperPairFlag());
        assertTrue(r2.getMateUnmappedFlag());
    }

    @Test
    public void txPrimaryWithRefXaAltAtSameLocusEmitsSplicedCigar()
    {
        // bwa primary on tx (spliced) with a ref XA at the same locus showing the alignment soft-clipped
        // around the junction. discriminator favours tx; XA emitted from the result should drop the ref
        // alt because it was Dropped by the discriminator.
        final SAMRecord r1 = newPrimary("read4", true, TX_CONTIG, 51, "100M");
        r1.setAttribute(XA_TAG, CHR_1 + ",+150,50M50S,0;");

        final SAMRecord r2 = newPrimary("read4", false, CHR_1, 700, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertEquals(CHR_1, r1.getReferenceName());
        assertEquals(150, r1.getAlignmentStart());
        assertEquals("50M100N50M", r1.getCigarString());
        // ref alt was Dropped — XA tag is null (only one kept alignment, which is the primary itself)
        assertNull(r1.getStringAttribute(XA_TAG));
    }

    // Issue #3: spliced records lifted from tx contigs should carry XS:A:+ / XS:A:- so downstream
    // RNA tools (Isofox) can interpret junctions strand-aware. Non-spliced records must NOT get XS:A.
    @Test
    public void splicedTxRecordsGetXsAStrand()
    {
        // forward-strand transcript (threeExonContig has strand=1): expect XS:A:+
        final SAMRecord r1 = newPrimary("readXsPlus", true, TX_CONTIG, 51, "100M");
        final SAMRecord r2 = newPrimary("readXsPlus", false, TX_CONTIG, 200, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        // r1 lifted to 50M100N50M → XS:A:+ on the spliced record
        assertTrue("expected N in lifted cigar", r1.getCigarString().contains("N"));
        assertEquals(Character.valueOf('+'), r1.getAttribute("XS"));
    }

    @Test
    public void nonSplicedRecordsDoNotGetXsA()
    {
        // ref-only alignment on chr1 (not from any tx contig), no N in cigar → no XS:A
        final SAMRecord r1 = newPrimary("readNoXs", true, CHR_1, 1000, "100M");
        final SAMRecord r2 = newPrimary("readNoXs", false, CHR_1, 1100, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertFalse("expected no N in lifted cigar", r1.getCigarString().contains("N"));
        assertNull(r1.getAttribute("XS"));
    }
}
