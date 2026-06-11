package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.SA_ATTRIBUTE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.redux.splice.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;
import com.hartwig.hmftools.redux.splice.rescue.RescueConfig;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

// Drives LiftBackGroupProcessor.processNameGroup with a capturing sink to assert emission/drop
// decisions. Focused on the orphaned-supplementary rule: a record that stays supplementary must carry
// an SA tag (REDUX dedup's FragmentCoords.fromRead reads the primary's coords from it), so a supp whose
// SA entries all fail to lift is dropped rather than emitted with a null SA.
public class LiftBackGroupProcessorTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    // exon spans on chr1: 100-199, 300-399, 500-549. introns: 200-299, 400-499. contig len 250.
    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    private static SAMRecord newPrimary(final String contig, final int pos, final String cigar)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("readX");
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setProperPairFlag(true);
        return record;
    }

    private static SAMRecord newSupplementary(final String contig, final int pos, final String cigar, final String saTag)
    {
        final SAMRecord record = newPrimary(contig, pos, cigar);
        record.setSupplementaryAlignmentFlag(true);
        record.setAttribute(SA_ATTRIBUTE, saTag);
        return record;
    }

    private static List<SAMRecord> process(final List<SAMRecord> group, final LiftBackStats stats)
    {
        return process(group, stats, null);
    }

    private static List<SAMRecord> process(
            final List<SAMRecord> group, final LiftBackStats stats, final JunctionRescueResolver rescueResolver)
    {
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        final LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, rescueResolver, null, null, null, null, 0, 0, stats);

        final List<SAMRecord> emitted = new ArrayList<>();
        processor.processNameGroup(group, new LiftedMateInfoCache(), (record, result) -> emitted.add(record));
        return emitted;
    }

    private static List<SAMRecord> processWithRef(final List<SAMRecord> group, final RefSequenceSource refSource)
    {
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        final LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, null, null, null, null, refSource, 0, 0, new LiftBackStats());

        final List<SAMRecord> emitted = new ArrayList<>();
        processor.processNameGroup(group, new LiftedMateInfoCache(), (record, result) -> emitted.add(record));
        return emitted;
    }

    @Test
    public void recomputesNmAgainstGenomicRefAndDropsMd()
    {
        // primary on the tx contig (exon1 100-199) lifts cleanly to chr1:100 50M. The genomic ref stub
        // returns all 'A'; the read has two trailing mismatches -> NM must be recomputed to 2, not carried
        // from the stale tx-contig NM:0, and MD must be dropped.
        final SAMRecord primary = newPrimary(TX_CONTIG, 1, "50M");
        final byte[] bases = new byte[50];
        Arrays.fill(bases, (byte) 'A');
        bases[48] = 'C';
        bases[49] = 'C';
        primary.setReadBases(bases);
        primary.setAttribute("NM", 0);
        primary.setAttribute("MD", "50");

        final RefSequenceSource refSource = (chromosome, posStart, posEnd) ->
        {
            final byte[] genomic = new byte[posEnd - posStart + 1];
            Arrays.fill(genomic, (byte) 'A');
            return genomic;
        };

        final List<SAMRecord> emitted = processWithRef(List.of(primary), refSource);

        assertEquals(1, emitted.size());
        final SAMRecord out = emitted.get(0);
        assertEquals(CHR_1, out.getReferenceName());
        assertEquals(Integer.valueOf(2), out.getIntegerAttribute("NM"));
        assertNull(out.getStringAttribute("MD"));
    }

    // rescue resolver with no annotated junctions: present (so the AS-unmap gate is active) but a no-op,
    // so a primary with no rescuable supps stays un-improved.
    private static JunctionRescueResolver noopRescue()
    {
        return new JunctionRescueResolver(Collections.emptySet(), RescueConfig.defaults());
    }

    @Test
    public void orphanedSupplementaryIsDropped()
    {
        // primary lifts cleanly; supp lifts too, but its only SA entry points at an out-of-range tx
        // position that fails to lift -> rewritten SA is null -> supp is dropped, not emitted.
        final SAMRecord primary = newPrimary(TX_CONTIG, 100, "50M");
        final SAMRecord supp = newSupplementary(TX_CONTIG, 110, "30M", TX_CONTIG + ",9999,+,30M,0,0;");

        final LiftBackStats stats = new LiftBackStats();
        final List<SAMRecord> emitted = process(List.of(primary, supp), stats);

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getSupplementaryAlignmentFlag());
        assertEquals(1, stats.orphanSuppsDropped());
    }

    @Test
    public void supplementaryWithLiftableSaIsKept()
    {
        // same shape, but the SA entry lifts -> supp is emitted with a rewritten genomic SA.
        final SAMRecord primary = newPrimary(TX_CONTIG, 100, "50M");
        final SAMRecord supp = newSupplementary(TX_CONTIG, 110, "30M", TX_CONTIG + ",100,+,50M,0,0;");

        final LiftBackStats stats = new LiftBackStats();
        final List<SAMRecord> emitted = process(List.of(primary, supp), stats);

        assertEquals(2, emitted.size());
        assertEquals(0, stats.orphanSuppsDropped());

        final SAMRecord emittedSupp = emitted.stream().filter(SAMRecord::getSupplementaryAlignmentFlag).findFirst().orElse(null);
        assertNotNull(emittedSupp);
        final String rewrittenSa = emittedSupp.getStringAttribute(SA_ATTRIBUTE);
        assertNotNull(rewrittenSa);
        assertTrue(rewrittenSa.startsWith(CHR_1 + ","));
    }

    @Test
    public void lowAsPrimaryUnmappedWhenLiftbackDidNotImprove()
    {
        // primary lifts but bwa scored it below the -T 30 floor and rescue (no junctions) can't improve it.
        final SAMRecord primary = newPrimary(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 20);

        final LiftBackStats stats = new LiftBackStats();
        final List<SAMRecord> emitted = process(List.of(primary), stats, noopRescue());

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
        assertEquals(1, stats.lowAsPrimariesUnmapped());
    }

    @Test
    public void highAsPrimaryKeptMapped()
    {
        final SAMRecord primary = newPrimary(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 60);

        final LiftBackStats stats = new LiftBackStats();
        final List<SAMRecord> emitted = process(List.of(primary), stats, noopRescue());

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
        assertEquals(0, stats.lowAsPrimariesUnmapped());
    }

    @Test
    public void lowAsPrimaryKeptWhenRescueDisabled()
    {
        // without a rescue resolver the AS-unmap gate is inactive, so a low-AS primary is left as-is.
        final SAMRecord primary = newPrimary(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 20);

        final LiftBackStats stats = new LiftBackStats();
        final List<SAMRecord> emitted = process(List.of(primary), stats);

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
        assertEquals(0, stats.lowAsPrimariesUnmapped());
    }
}
