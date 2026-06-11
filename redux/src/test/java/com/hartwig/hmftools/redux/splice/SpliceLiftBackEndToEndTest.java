package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.splice.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.redux.splice.LiftBackRecordOps.applyResultToRecord;
import static com.hartwig.hmftools.redux.splice.LiftBackRecordOps.toLiftedMateInfo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

// Black-box end-to-end tests for the full SpliceLiftBack pipeline (resolver + apply + SA rewrite +
// mate patch). Tests are written in terms of bwa input and expected BAM output — no inner-class access.
public class SpliceLiftBackEndToEndTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;
    private static final String XA_TAG = "XA";

    // exon spans: 100-199, 300-399, 500-549; introns: 200-299, 400-499; contig len 250.
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

    private static List<SAMRecord> processFragment(final List<ContigEntry> entries, final List<SAMRecord> records)
    {
        final LiftBackResolver resolver = new LiftBackResolver(entries);
        final LiftedMateInfoCache cache = new LiftedMateInfoCache();

        for(final SAMRecord record : records)
        {
            if(!record.getReadPairedFlag() || record.isSecondaryOrSupplementary())
                continue;
            final LiftBackResult result = resolver.resolve(record);
            cache.recordPrimaryAlignment(record.getReadName(), record.getFirstOfPairFlag(), toLiftedMateInfo(record, result));
        }

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
        // R1 spans exon1->exon2; R2 starts at tx 197 so its leading exon2 anchor is 4M (above
        // trimMicroAnchors threshold of 3 — tx 200 would produce a 1M head that would be dropped).
        final SAMRecord r1 = newPrimary("read1", true, TX_CONTIG, 51, "100M");
        final SAMRecord r2 = newPrimary("read1", false, TX_CONTIG, 197, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertEquals(CHR_1, r1.getReferenceName());
        assertEquals(150, r1.getAlignmentStart());
        assertEquals("50M100N50M", r1.getCigarString());
        assertEquals(CHR_1, r1.getMateReferenceName());
        assertEquals(396, r1.getMateAlignmentStart());

        assertEquals(CHR_1, r2.getReferenceName());
        assertEquals(396, r2.getAlignmentStart());
        assertTrue(r2.getCigarString().contains("N"));
        assertEquals(150, r2.getMateAlignmentStart());
    }

    @Test
    public void tinyJunctionAnchorKeptNotClampedToSoftclip()
    {
        // R2 starts at tx 200 (last base of exon2), producing a 1bp exon2 anchor. With
        // ANNOTATED_JUNCTION_MIN_ANCHOR_BP=1 this is kept as a real junction (1M...N...M) matching STAR,
        // not rolled into a softclip. Guards against silently re-raising the floor to 3.
        final SAMRecord r1 = newPrimary("read3", true, TX_CONTIG, 51, "100M");
        final SAMRecord r2 = newPrimary("read3", false, TX_CONTIG, 200, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        // R2: 1bp exon2 anchor + 100N intron2 + 49bp exon3.
        assertEquals(CHR_1, r2.getReferenceName());
        assertEquals(399, r2.getAlignmentStart());
        assertEquals("1M100N49M", r2.getCigarString());
    }

    @Test
    public void supplementarySaTagRewrittenToGenomicCoords()
    {
        // Supp starts at tx 197 so its exon2 anchor is 4M (above trimMicroAnchors threshold of 3).
        final SAMRecord r1 = newPrimary("read2", true, TX_CONTIG, 51, "100M");
        r1.setAttribute(SA_ATTRIBUTE, TX_CONTIG + ",197,+,50M,60,0;");

        final SAMRecord r1Supp = newSupplementary("read2", true, TX_CONTIG, 197, "50M");
        r1Supp.setAttribute(SA_ATTRIBUTE, TX_CONTIG + ",51,+,100M,60,0;");

        final SAMRecord r2 = newPrimary("read2", false, CHR_1, 600, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r1Supp, r2));

        final String r1Sa = r1.getStringAttribute(SA_ATTRIBUTE);
        assertTrue("expected lifted SA entry on chr1, got: " + r1Sa, r1Sa != null && r1Sa.startsWith(CHR_1 + ","));
        assertFalse(r1Sa.contains(TX_CONTIG));

        assertEquals(CHR_1, r1Supp.getReferenceName());
        assertEquals(396, r1Supp.getAlignmentStart());

        assertEquals(CHR_1, r1.getMateReferenceName());
        assertEquals(600, r1.getMateAlignmentStart());
        assertEquals(CHR_1, r1Supp.getMateReferenceName());
        assertEquals(600, r1Supp.getMateAlignmentStart());
    }

    @Test
    public void unliftablePrimaryFlagsRecordUnmappedAndClearsMatesProperPair()
    {
        // R1 starts past altEnd — translation fails, emitted as unmapped.
        final SAMRecord r1 = newPrimary("read3", true, TX_CONTIG, 251, "50M");
        final SAMRecord r2 = newPrimary("read3", false, CHR_1, 600, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertTrue(r1.getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, r1.getCigarString());
        assertNull(r1.getStringAttribute(XA_TAG));

        assertEquals(CHR_1, r2.getReferenceName());
        assertFalse(r2.getProperPairFlag());
        assertTrue(r2.getMateUnmappedFlag());
    }

    @Test
    public void txPrimaryWithRefXaAltAtSameLocusEmitsSplicedCigar()
    {
        // bwa primary on tx with a softclipped ref XA at the same locus; discriminator favours tx and drops the ref alt.
        final SAMRecord r1 = newPrimary("read4", true, TX_CONTIG, 51, "100M");
        r1.setAttribute(XA_TAG, CHR_1 + ",+150,50M50S,0;");

        final SAMRecord r2 = newPrimary("read4", false, CHR_1, 700, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertEquals(CHR_1, r1.getReferenceName());
        assertEquals(150, r1.getAlignmentStart());
        assertEquals("50M100N50M", r1.getCigarString());
        assertNull(r1.getStringAttribute(XA_TAG));
    }

    // Spliced tx records must carry XS:A:+/- for strand-aware junction interpretation (e.g. Isofox).
    // Non-spliced records must NOT get XS:A.
    @Test
    public void splicedTxRecordsGetXsAStrand()
    {
        final SAMRecord r1 = newPrimary("readXsPlus", true, TX_CONTIG, 51, "100M");
        final SAMRecord r2 = newPrimary("readXsPlus", false, TX_CONTIG, 200, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertTrue("expected N in lifted cigar", r1.getCigarString().contains("N"));
        assertEquals(Character.valueOf('+'), r1.getAttribute("XS"));
    }

    @Test
    public void nonSplicedRecordsDoNotGetXsA()
    {
        final SAMRecord r1 = newPrimary("readNoXs", true, CHR_1, 1000, "100M");
        final SAMRecord r2 = newPrimary("readNoXs", false, CHR_1, 1100, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertFalse("expected no N in lifted cigar", r1.getCigarString().contains("N"));
        assertNull(r1.getAttribute("XS"));
    }
}
