package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.refSource;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.secondMateRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.tars.liftback.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.RescueConfig;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// Drives LiftBackGroupProcessor.processNameGroup with a capturing sink to assert emission/drop
// decisions. Focused on the orphaned-supplementary rule: a record that stays supplementary must carry
// an SA tag (REDUX dedup's FragmentCoords.fromRead reads the primary's coords from it), so a supp whose
// SA entries all fail to lift is dropped rather than emitted with a null SA.
public class LiftBackGroupProcessorTest
{
    private static List<SAMRecord> process(final List<SAMRecord> group, final LiftBackStats stats)
    {
        return process(group, stats, null, null, 0, 0);
    }

    private static List<SAMRecord> process(
            final List<SAMRecord> group, final LiftBackStats stats, final JunctionRescueResolver rescueResolver)
    {
        return process(group, stats, rescueResolver, null, 0, 0);
    }

    // Single harness for all LiftBackGroupProcessor scenarios: the only thing tests vary is which ctor arg is
    // non-default, so a ctor change touches one call site.
    private static List<SAMRecord> process(
            final List<SAMRecord> group, final LiftBackStats stats, final JunctionRescueResolver rescueResolver,
            final RefSequenceSource refSource, final int unmapAboveNh, final int unmapBelowMapq)
    {
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        final LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, rescueResolver, null, null, null, refSource, unmapAboveNh, unmapBelowMapq, stats);

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
        final SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        primary.setReadBases(bases("A".repeat(48) + "CC"));   // two trailing mismatches vs the all-'A' genomic stub
        primary.setAttribute("NM", 0);
        primary.setAttribute("MD", "50");

        final String chr1Bases = "A".repeat(600);
        final List<SAMRecord> emitted = process(List.of(primary), new LiftBackStats(), null, refSource(CHR_1, chr1Bases), 0, 0);

        assertEquals(1, emitted.size());
        final SAMRecord out = emitted.get(0);
        assertEquals(CHR_1, out.getReferenceName());
        assertEquals(Integer.valueOf(2), out.getIntegerAttribute("NM"));
        assertNull(out.getStringAttribute("MD"));
    }

    @Test
    public void emittedPrimaryTaggedWithLocusCountNh()
    {
        // single-locus ref-only primary -> NH = 1.
        final SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");

        final List<SAMRecord> emitted = process(List.of(primary), new LiftBackStats());

        assertEquals(1, emitted.size());
        assertEquals(Integer.valueOf(1), emitted.get(0).getIntegerAttribute("NH"));
    }

    @Test
    public void unmapBelowMapqUnmapsLowMapqPrimary()
    {
        final SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        primary.setMappingQuality(10);

        final List<SAMRecord> emitted = process(List.of(primary), new LiftBackStats(), null, null, 0, 20);

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void mapqAtOrAboveFloorKeptMapped()
    {
        final SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        primary.setMappingQuality(20);

        final List<SAMRecord> emitted = process(List.of(primary), new LiftBackStats(), null, null, 0, 20);

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void pairedMatesPatchedAgainstPerGroupCache()
    {
        // both mates of one fragment in the group: /1 lifts to chr1:100, /2 (exon2) to chr1:300. The
        // per-group cache must let each mate's fields point at the other's lifted coords -- this is the
        // single-pass per-group correctness property that replaced the whole-sample pass-1 cache.
        final SAMRecord mate1 = primaryRecord(TX_CONTIG, 1, "50M");
        final SAMRecord mate2 = secondMateRecord(TX_CONTIG, 101, "50M");

        final List<SAMRecord> emitted = process(List.of(mate1, mate2), new LiftBackStats());

        assertEquals(2, emitted.size());
        final SAMRecord out1 = emitted.stream().filter(SAMRecord::getFirstOfPairFlag).findFirst().orElseThrow();
        final SAMRecord out2 = emitted.stream().filter(record -> !record.getFirstOfPairFlag()).findFirst().orElseThrow();

        assertEquals(CHR_1, out1.getReferenceName());
        assertEquals(100, out1.getAlignmentStart());
        assertEquals(CHR_1, out2.getReferenceName());
        assertEquals(300, out2.getAlignmentStart());

        assertEquals(CHR_1, out1.getMateReferenceName());
        assertEquals(300, out1.getMateAlignmentStart());
        assertEquals(CHR_1, out2.getMateReferenceName());
        assertEquals(100, out2.getMateAlignmentStart());
    }

    @Test
    public void duplicateSupplementariesCollapsed()
    {
        // two supps lifting to the same (chrom, pos, cigar, strand) collapse to one -- bwa can emit the same
        // junction across multiple tx contigs.
        final SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        final SAMRecord supp1 = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",1,+,50M,0,0;");
        final SAMRecord supp2 = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",1,+,50M,0,0;");

        final List<SAMRecord> emitted = process(List.of(primary, supp1, supp2), new LiftBackStats());

        assertEquals(2, emitted.size());
        assertEquals(1, emitted.stream().filter(SAMRecord::getSupplementaryAlignmentFlag).count());
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
        final SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        final SAMRecord supp = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",9999,+,30M,0,0;");

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
        final SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        final SAMRecord supp = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",100,+,50M,0,0;");

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
        final SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
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
        final SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
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
        final SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 20);

        final LiftBackStats stats = new LiftBackStats();
        final List<SAMRecord> emitted = process(List.of(primary), stats);

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
        assertEquals(0, stats.lowAsPrimariesUnmapped());
    }
}
