package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.pairedUnmappedRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.refGenome;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.secondMateRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.unpairedPrimaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.unpairedSupplementaryRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;
import com.hartwig.hmftools.tars.liftback.features.GenomicAlignmentScorer;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryConfig;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// Drives LiftBackGroupProcessor.processNameGroup with a capturing consumer to assert which records are emitted,
// dropped or unmapped. The orphaned-supplementary rule is the load-bearing one: a record that stays supplementary
// must carry an SA tag, because REDUX dedup reads the primary's coords from it, so a supp whose SA entries all fail
// to lift is dropped rather than emitted with a null SA.
public class LiftBackGroupProcessorTest
{
    // Every scenario varies only which collaborator is non-null, so one full harness keeps a ctor change to a
    // single call site.
    private static List<SAMRecord> process(
            final List<SAMRecord> group, final SupplementaryResolver supplementaryResolver,
            final OverhangGate overhangGate, final RefGenomeInterface refGenome, final ExcludedRegions excludedRegions)
    {
        LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                new LiftBackDiscriminator(List.of(threeExonContig())),
                supplementaryResolver, overhangGate, new GenomicAlignmentScorer(refGenome), refGenome, excludedRegions);

        List<SAMRecord> emitted = new ArrayList<>();
        processor.processNameGroup(group, emitted::add);
        return emitted;
    }

    private static List<SAMRecord> process(final List<SAMRecord> group)
    {
        return process(group, null, null, null, null);
    }

    private static List<SAMRecord> process(final List<SAMRecord> group, final SupplementaryResolver supplementaryResolver)
    {
        return process(group, supplementaryResolver, null, null, null);
    }

    private static ExcludedRegions excludedRegion(final String chromosome, final int start, final int end)
    {
        // mutable: the ExcludedRegions ctor sorts the list in place
        Map<String, List<BaseRegion>> regions = new HashMap<>();
        regions.put(chromosome, new ArrayList<>(List.of(new BaseRegion(start, end))));
        return new ExcludedRegions(regions);
    }

    // Supplementary resolver with no annotated junctions: present, so the AS-unmap gate is active, but a no-op,
    // so a primary with no resolvable supps stays un-improved.
    private static SupplementaryResolver noopSupplementary()
    {
        return new SupplementaryResolver(Collections.emptySet(), SupplementaryConfig.defaults());
    }

    @Test
    public void testPrimaryLiftingIntoExcludedRegionIsUnmapped()
    {
        // tx primary (exon1) lifts to chr1:100; an excluded region covering it unmaps the read REDUX-style:
        // kept in the output but flagged unmapped with no cigar, not aligned in the excluded zone.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");

        List<SAMRecord> emitted = process(List.of(primary), null, null, null, excludedRegion(CHR_1, 50, 300));

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, emitted.get(0).getCigarString());
    }

    @Test
    public void testPrimaryOutsideExcludedRegionIsKept()
    {
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");   // chr1:100

        List<SAMRecord> emitted = process(List.of(primary), null, null, null, excludedRegion(CHR_1, 5000, 6000));

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void testRecomputesNmAgainstGenomicRefAndDropsMd()
    {
        // primary on the tx contig (exon1 100-199) lifts cleanly to chr1:100 50M. The genomic ref stub returns all
        // 'A'; the read has two trailing mismatches, so NM must be recomputed to 2 rather than carried from the
        // stale tx-contig NM:0, and MD must be dropped.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        primary.setReadBases(bases("A".repeat(48) + "CC"));
        primary.setAttribute("NM", 0);
        primary.setAttribute("MD", "50");

        List<SAMRecord> emitted = process(List.of(primary), null, null, refGenome(CHR_1, "A".repeat(600)), null);

        assertEquals(1, emitted.size());
        SAMRecord out = emitted.get(0);
        assertEquals(CHR_1, out.getReferenceName());
        assertEquals(Integer.valueOf(2), out.getIntegerAttribute("NM"));
        assertNull(out.getStringAttribute("MD"));
    }

    @Test
    public void testEmittedPrimaryTaggedWithLocusCountNh()
    {
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");   // single-locus ref-only primary

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        assertEquals(Integer.valueOf(1), emitted.get(0).getIntegerAttribute("NH"));
    }

    @Test
    public void testOverCapGenomicPrimaryMapQuality0NoXaUnmapped()
    {
        // a GENOMIC primary emitted MAPQ 0 with no XA maps past the XA cap (75+ distinct genomic loci), so it is
        // unmapped even though, with no XA, the discriminator sees a single locus and would otherwise bump to 60.
        SAMRecord primary = primaryRecord(CHR_1, 100, "50M");
        primary.setMappingQuality(0);

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void testOverCapTxPrimaryMapQuality0NoXaKept()
    {
        // a TX-CONTIG primary with MAPQ 0 and no XA hit 75+ transcript contigs of one gene, which all lift to one
        // genomic locus - the over-cap rule is ref-only-gated and must not unmap it.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        primary.setMappingQuality(0);

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
        assertEquals(CHR_1, emitted.get(0).getReferenceName());
    }

    @Test
    public void testMapQuality0WithXaKeptMapped()
    {
        // MAPQ 0 but XA present is an ordinary few-way multimapper within the cap; keep and lift it.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        primary.setMappingQuality(0);
        primary.setAttribute("XA", CHR_1 + ",+5000,50M,0;");

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void testPairedMatesPatchedAgainstPerGroupPair()
    {
        // both mates of one fragment in the group: /1 lifts to chr1:100, /2 (exon2) to chr1:300. The per-group
        // pair must let each mate's fields point at the other's lifted coords - the single-pass correctness
        // property that replaced the whole-sample first pass.
        SAMRecord mate1 = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord mate2 = secondMateRecord(TX_CONTIG, 101, "50M");

        List<SAMRecord> emitted = process(List.of(mate1, mate2));

        assertEquals(2, emitted.size());
        SAMRecord out1 = emitted.stream().filter(SAMRecord::getFirstOfPairFlag).findFirst().orElseThrow();
        SAMRecord out2 = emitted.stream().filter(record -> !record.getFirstOfPairFlag()).findFirst().orElseThrow();

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
    public void testOneUnmappedMateIsPlacedBesideMappedMateReduxStyle()
    {
        // /1 is deliberately unmapped by TARS's genomic over-cap rule while /2 remains mapped. REDUX represents
        // that fragment by parking /1 at /2's genomic coordinate; /2 points its unmapped-mate fields back to itself.
        SAMRecord mate1 = primaryRecord(CHR_1, 100, "50M");
        mate1.setMappingQuality(0);
        SAMRecord mate2 = secondMateRecord(CHR_1, 600, "50M");

        List<SAMRecord> emitted = process(List.of(mate1, mate2));

        assertEquals(2, emitted.size());
        SAMRecord out1 = emitted.stream().filter(SAMRecord::getFirstOfPairFlag).findFirst().orElseThrow();
        SAMRecord out2 = emitted.stream().filter(record -> !record.getFirstOfPairFlag()).findFirst().orElseThrow();

        assertTrue(out1.getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, out1.getCigarString());
        assertEquals(0, out1.getMappingQuality());
        assertEquals("1:100", out1.getStringAttribute("UM"));
        assertEquals(out2.getReferenceName(), out1.getReferenceName());
        assertEquals(out2.getAlignmentStart(), out1.getAlignmentStart());
        assertFalse(out1.getMateUnmappedFlag());
        assertEquals(out2.getReferenceName(), out1.getMateReferenceName());
        assertEquals(out2.getAlignmentStart(), out1.getMateAlignmentStart());

        assertFalse(out2.getReadUnmappedFlag());
        assertTrue(out2.getMateUnmappedFlag());
        assertEquals(out2.getReferenceName(), out2.getMateReferenceName());
        assertEquals(out2.getAlignmentStart(), out2.getMateAlignmentStart());
        assertFalse(out1.getProperPairFlag());
        assertFalse(out2.getProperPairFlag());
        assertEquals(0, out1.getInferredInsertSize());
        assertEquals(0, out2.getInferredInsertSize());
    }

    @Test
    public void testSingleEndPrimaryLiftedWithoutMatePatching()
    {
        // Ultima single-end: an unpaired primary must lift with no mate patching, and without the
        // getFirstOfPairFlag() throw the per-group mate refresh used to hit on unpaired input.
        SAMRecord primary = unpairedPrimaryRecord(TX_CONTIG, 1, "50M");

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        SAMRecord out = emitted.get(0);
        assertFalse(out.getReadPairedFlag());
        assertFalse(out.getReadUnmappedFlag());
        assertEquals(CHR_1, out.getReferenceName());
        assertEquals(100, out.getAlignmentStart());
        assertEquals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, out.getMateReferenceName());
        assertEquals(0, out.getInferredInsertSize());
    }

    @Test
    public void testSingleEndPrimaryWithSupplementaryLifted()
    {
        // unpaired primary plus its supplementary both lift and emit; exercises the mate refresh with a supp in
        // the group, which is where the unpaired getFirstOfPairFlag() call sat.
        SAMRecord primary = unpairedPrimaryRecord(TX_CONTIG, 100, "50M");
        SAMRecord supp = unpairedSupplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",100,+,50M,0,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp));

        assertEquals(2, emitted.size());
        assertTrue(emitted.stream().noneMatch(SAMRecord::getReadPairedFlag));
        assertEquals(1, emitted.stream().filter(SAMRecord::getSupplementaryAlignmentFlag).count());
    }

    @Test
    public void testDuplicateSupplementariesCollapsed()
    {
        // two supps lifting to the same chromosome, position, cigar and strand collapse to one - bwa can emit the
        // same junction across multiple tx contigs.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord supp1 = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",1,+,50M,0,0;");
        SAMRecord supp2 = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",1,+,50M,0,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp1, supp2));

        assertEquals(2, emitted.size());
        assertEquals(1, emitted.stream().filter(SAMRecord::getSupplementaryAlignmentFlag).count());
    }

    @Test
    public void testOrphanedSupplementaryIsDropped()
    {
        // primary lifts cleanly; the supp lifts too, but its only SA entry points at an out-of-range tx position
        // that fails to lift, so the rewritten SA is null and the supp is dropped rather than emitted.
        SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        SAMRecord supp = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",9999,+,30M,0,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp));

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getSupplementaryAlignmentFlag());
    }

    @Test
    public void testSupplementaryWithLiftableSaIsKept()
    {
        // same shape, but the SA entry lifts, so the supp is emitted with a rewritten genomic SA.
        SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        SAMRecord supp = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",100,+,50M,0,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp));

        assertEquals(2, emitted.size());

        SAMRecord emittedSupp = emitted.stream().filter(SAMRecord::getSupplementaryAlignmentFlag).findFirst().orElse(null);
        assertNotNull(emittedSupp);
        String rewrittenSa = emittedSupp.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE);
        assertNotNull(rewrittenSa);
        assertTrue(rewrittenSa.startsWith(CHR_1 + ","));
    }

    @Test
    public void testLowAsPrimaryUnmappedWhenLiftbackDidNotImprove()
    {
        // primary lifts but bwa scored it below the -T 30 floor, and supplementary resolve (no junctions) cannot
        // improve it.
        SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 20);

        List<SAMRecord> emitted = process(List.of(primary), noopSupplementary());

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void testHighAsPrimaryKeptMapped()
    {
        SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 60);

        List<SAMRecord> emitted = process(List.of(primary), noopSupplementary());

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void testLowAsPrimaryKeptWithoutSupplementaryResolver()
    {
        // with no supplementary resolver the AS-unmap gate is inactive, so a low-AS primary is left as-is.
        SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 20);

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void testLowAsPrimaryUnmapDropsItsSupplementary()
    {
        // The AS-floor unmap happens during decide, so the supp-orphan gate observes the unmapped primary and drops
        // the supplementary. Firing it only at emit would leave the chosen primary mapped, so the supp would
        // survive with an SA referencing a primary emitted as unmapped - a dangling locus for REDUX FragmentCoords.
        // The supp's SA lifts cleanly, so without the decide-time unmap this would emit two records, not one.
        SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        primary.setAttribute("AS", 20);
        SAMRecord supp = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",100,+,50M,0,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp), noopSupplementary());

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
        assertFalse(emitted.get(0).getSupplementaryAlignmentFlag());
    }

    @Test
    public void testMultimapperMatesColocatedAtOneLocus()
    {
        // Real-shape multimapper (the TERT / chr10 field case): a fully-overlapping short-insert pair whose two
        // mates map equally well (NM 0, tied genomic score) to two genomic copies. Through the full group processor
        // the fragment must not be torn across the copies - both mates land on one chromosome, the other in XA.
        String seq = "ACGTTGCA".repeat(7).substring(0, 55);   // both copies carry this exact sequence
        RefGenomeInterface ref = new TestGenome()
                .with("chr5", 300, 'A').with("chr10", 300, 'A')
                .set("chr5", 100, seq).set("chr10", 100, seq)
                .asRefGenome();

        SAMRecord mate1 = primaryRecord("frag", "chr5", 100, "55M");   // read /1, reverse strand
        mate1.setMappingQuality(0);
        mate1.setReadNegativeStrandFlag(true);
        mate1.setProperPairFlag(false);
        mate1.setReadBases(bases(seq));
        mate1.setAttribute("XA", "chr10,-100,55M,0;");

        SAMRecord mate2 = secondMateRecord("frag", "chr10", 100, "55M");   // read /2, forward strand
        mate2.setMappingQuality(0);
        mate2.setProperPairFlag(false);
        mate2.setReadBases(bases(seq));
        mate2.setAttribute("XA", "chr5,+100,55M,0;");

        List<SAMRecord> emitted = process(List.of(mate1, mate2), null, new OverhangGate(ref), ref, null);

        assertEquals(2, emitted.size());
        SAMRecord out1 = emitted.stream().filter(SAMRecord::getFirstOfPairFlag).findFirst().orElseThrow();
        SAMRecord out2 = emitted.stream().filter(record -> !record.getFirstOfPairFlag()).findFirst().orElseThrow();

        assertFalse(out1.getReadUnmappedFlag());
        assertFalse(out2.getReadUnmappedFlag());
        assertEquals("both mates land on one chromosome", out1.getReferenceName(), out2.getReferenceName());
        assertEquals(out1.getReferenceName(), out1.getMateReferenceName());
        assertEquals(out2.getReferenceName(), out2.getMateReferenceName());

        String otherCopy = out1.getReferenceName().equals("chr5") ? "chr10" : "chr5";
        assertTrue(
                "other copy carried in XA",
                out1.getStringAttribute("XA") != null && out1.getStringAttribute("XA").contains(otherCopy));
    }

    @Test
    public void testFirstMateUsesSecondMateCandidates()
    {
        String sequence = "ACGT".repeat(13);
        RefGenomeInterface ref = new TestGenome()
                .with("chr5", 500, 'A').with("chr10", 500, 'A')
                .set("chr5", 100, sequence).set("chr10", 100, sequence)
                .asRefGenome();

        SAMRecord first = primaryRecord("paired", "chr5", 100, "52M");
        first.setMappingQuality(0);
        first.setReadBases(bases(sequence));
        first.setAttribute("XA", "chr10,+100,52M,0;");

        SAMRecord second = secondMateRecord("paired", "chr10", 120, "52M");
        second.setReadBases(bases(sequence));

        List<SAMRecord> emitted = process(List.of(first, second), null, new OverhangGate(ref), ref, null);

        assertEquals("chr10", emitted.stream().filter(SAMRecord::getFirstOfPairFlag).findFirst().orElseThrow().getReferenceName());
    }

    @Test
    public void testProcessorMergesSupplementaryChainBeforeDiscrimination()
    {
        SAMRecord primary = primaryRecord(CHR_1, 1000, "50M101S");
        SAMRecord middle = supplementaryRecord(CHR_1, 2000, "50S60M41S", CHR_1 + ",1000,+,50M101S,0,0;");
        SAMRecord last = supplementaryRecord(CHR_1, 3000, "110S41M", CHR_1 + ",1000,+,50M101S,0,0;");
        primary.setReadBases(bases("A".repeat(151)));
        middle.setReadBases(bases("A".repeat(151)));
        last.setReadBases(bases("A".repeat(151)));
        SupplementaryResolver resolver = new SupplementaryResolver(
                Set.of(new ChrBaseRegion(CHR_1, 1050, 1999), new ChrBaseRegion(CHR_1, 2060, 2999)),
                SupplementaryConfig.defaults());

        List<SAMRecord> emitted = process(List.of(primary, middle, last), resolver);

        assertEquals(1, emitted.size());
        assertEquals("50M950N60M940N41M", emitted.get(0).getCigarString());
    }

    @Test
    public void testGenomicSupplementaryMergeGetsXsFromAnnotation()
    {
        SAMRecord primary = primaryRecord(CHR_1, 150, "50M101S");
        SAMRecord supplementary = supplementaryRecord(CHR_1, 300, "50S101M", CHR_1 + ",150,+,50M101S,0,0;");
        primary.setReadBases(bases("A".repeat(151)));
        supplementary.setReadBases(bases("A".repeat(151)));

        List<SAMRecord> emitted = process(List.of(primary, supplementary), contigSupplementary());

        assertEquals(1, emitted.size());
        assertEquals("50M100N101M", emitted.get(0).getCigarString());
        assertEquals(Character.valueOf('+'), emitted.get(0).getAttribute("XS"));
    }

    @Test
    public void testFinalGenomicScoreControlsAlignmentScoreFloor()
    {
        RefGenomeInterface ref = new TestGenome().with(CHR_1, 500, 'A').asRefGenome();
        SAMRecord primary = primaryRecord(CHR_1, 1, "20M100N3M48S");
        primary.setReadBases(bases("C".repeat(71)));
        primary.setAttribute("AS", 60);

        List<SAMRecord> emitted = process(
                List.of(primary), noopSupplementary(), new OverhangGate(ref), ref, null);

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
    }

    // bwa emits a pair with both ends unmapped as flags 77/141 - paired, read unmapped, mate unmapped. The mate-unmapped
    // bit must survive the lift: clearing it leaves a record claiming a mapped mate with RNEXT unset, which REDUX rejects
    // with htsjdk's INVALID_FLAG_MATE_UNMAPPED.
    @Test
    public void testBothEndsUnmappedPairKeepsMateUnmappedFlag()
    {
        SAMRecord first = pairedUnmappedRecord("frag1", true);
        first.setMateUnmappedFlag(true);

        SAMRecord second = pairedUnmappedRecord("frag1", false);
        second.setMateUnmappedFlag(true);

        List<SAMRecord> emitted = process(List.of(first, second));

        assertEquals(2, emitted.size());
        for(SAMRecord record : emitted)
        {
            String side = record.getFirstOfPairFlag() ? "first" : "second";
            assertEquals(side + " RNEXT", NO_CHROMOSOME_NAME, record.getMateReferenceName());
            assertEquals(side + " flags", record.getFirstOfPairFlag() ? 77 : 141, record.getFlags());
            assertTrue(side + " read should stay unmapped", record.getReadUnmappedFlag());
            assertTrue(side + " mate-unmapped flag should survive the lift", record.getMateUnmappedFlag());
        }
    }

    // Isolates the mate patch from the group processor: an unrecorded mate must leave the mate-unmapped bit set.
    @Test
    public void testPatchMateFieldsKeepsMateUnmappedForUnplacedMate()
    {
        SAMRecord record = pairedUnmappedRecord("frag1", true);
        record.setMateUnmappedFlag(true);

        new LiftedMatePair().patchMateFields(record);

        assertEquals("flags", 77, record.getFlags());
    }

    // Junctions taken from the same sidecar entry the discriminator lifts against, so the intron coords and the
    // chromosome key match what the lift emits: chr1 introns 200-299 and 400-499 between the three exons.
    private static SupplementaryResolver contigSupplementary()
    {
        return new SupplementaryResolver(
                EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())), null, SupplementaryConfig.defaults());
    }

}
