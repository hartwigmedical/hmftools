package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.applyResultToRecord;
import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.toLiftedMateInfo;
import static com.hartwig.hmftools.tars.liftback.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.assertLifted;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.secondMateRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// Black-box end-to-end tests for the full SpliceLiftBack pipeline (resolver + apply + SA rewrite +
// mate patch). Tests are written in terms of bwa input and expected BAM output — no inner-class access.
public class SpliceLiftBackEndToEndTest
{
    private static final String XA_TAG = "XA";

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
    public void tinyJunctionAnchorKeptNotClampedToSoftclip()
    {
        // R2 starts at tx 200 (last base of exon2), producing a 1bp exon2 anchor. With
        // ANNOTATED_JUNCTION_MIN_ANCHOR_BP=1 this is kept as a real junction (1M...N...M),
        // not rolled into a softclip. Guards against silently re-raising the floor to 3.
        final SAMRecord r1 = primaryRecord("read3", TX_CONTIG, 51, "100M");
        final SAMRecord r2 = secondMateRecord("read3", TX_CONTIG, 200, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        // R2: 1bp exon2 anchor + 100N intron2 + 49bp exon3.
        assertLifted(r2, CHR_1, 399, "1M100N49M");
    }

    @Test
    public void txPrimaryWithRefXaAltAtSameLocusEmitsSplicedCigar()
    {
        // bwa primary on tx with a softclipped ref XA at the same locus; discriminator favours tx and drops the ref alt.
        final SAMRecord r1 = primaryRecord("read4", TX_CONTIG, 51, "100M");
        r1.setAttribute(XA_TAG, CHR_1 + ",+150,50M50S,0;");

        final SAMRecord r2 = secondMateRecord("read4", CHR_1, 700, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertLifted(r1, CHR_1, 150, "50M100N50M");
        assertNull(r1.getStringAttribute(XA_TAG));
    }

    // Spliced tx records must carry XS:A:+/- for strand-aware junction interpretation (e.g. Isofox).
    // Non-spliced records must NOT get XS:A.
    @Test
    public void splicedTxRecordsGetXsAStrand()
    {
        final SAMRecord r1 = primaryRecord("readXsPlus", TX_CONTIG, 51, "100M");
        final SAMRecord r2 = secondMateRecord("readXsPlus", TX_CONTIG, 200, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertTrue("expected N in lifted cigar", r1.getCigarString().contains("N"));
        assertEquals(Character.valueOf('+'), r1.getAttribute("XS"));
    }

    @Test
    public void nonSplicedRecordsDoNotGetXsA()
    {
        final SAMRecord r1 = primaryRecord("readNoXs", CHR_1, 1000, "100M");
        final SAMRecord r2 = secondMateRecord("readNoXs", CHR_1, 1100, "50M");

        processFragment(List.of(threeExonContig()), List.of(r1, r2));

        assertFalse("expected no N in lifted cigar", r1.getCigarString().contains("N"));
        assertNull(r1.getAttribute("XS"));
    }
}
