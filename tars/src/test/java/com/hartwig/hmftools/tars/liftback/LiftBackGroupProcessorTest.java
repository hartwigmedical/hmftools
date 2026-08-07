package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
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
import com.hartwig.hmftools.tars.liftback.features.SoftClipExtender;
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
                supplementaryResolver, overhangGate, new SoftClipExtender(refGenome), refGenome, excludedRegions);

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

    // Junctions taken from the same sidecar entry the discriminator lifts against, so the intron coords and the
    // chromosome key match what the lift emits: chr1 introns 200-299 and 400-499 between the three exons.
    private static SupplementaryResolver contigSupplementary()
    {
        return new SupplementaryResolver(
                EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())), null, SupplementaryConfig.defaults());
    }

    @Test
    public void testBoundarySnapReachesBothRecordsOfASplitRead()
    {
        // The fusion shape: 51M ends at 200, one base into the annotated intron, and the supp starts at 499, one base
        // before the exon at 500. Their read coverage leaves a gap so no merge can carry either correction. Asserted on
        // the emitted records rather than on the resolver, because the snap has to survive the overhang gate to be worth
        // having, and on both records because a fusion's two junction ends come from two of them.
        SAMRecord primary = primaryRecord(CHR_1, 150, "51M100S");
        SAMRecord supp = supplementaryRecord(CHR_1, 499, "100S51M", CHR_1 + ",150,+,51M100S,60,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp), contigSupplementary());

        assertEquals(2, emitted.size());
        assertEquals("50M101S", emitted.get(0).getCigarString());
        assertEquals(150, emitted.get(0).getAlignmentStart());
        assertEquals("101S50M", emitted.get(1).getCigarString());
        assertEquals("left retraction moves the start", 500, emitted.get(1).getAlignmentStart());
    }

    @Test
    public void testBoundarySnapDoesNotRescueALowAlignmentScorePrimary()
    {
        // AS 20 is under the floor, so the primary unmaps and its supp is dropped as an orphan. The snap must not change
        // that: it shortens the alignment, so bwa's stale AS is if anything generous. Suppression is for the merge and
        // collapse passes, which lengthen it.
        SAMRecord primary = primaryRecord(CHR_1, 150, "51M100S");
        primary.setAttribute("AS", 20);
        SAMRecord supp = supplementaryRecord(CHR_1, 499, "100S51M", CHR_1 + ",150,+,51M100S,60,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp), contigSupplementary());

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
        assertFalse(emitted.get(0).getSupplementaryAlignmentFlag());
    }

    @Test
    public void testNoBoundarySnapForAReadWithNoSupplementary()
    {
        // Same over-extended boundary as above, but with no partner record: this is as likely to be intron retention as
        // over-extension, so the bases stay aligned.
        SAMRecord primary = primaryRecord(CHR_1, 150, "51M100S");

        List<SAMRecord> emitted = process(List.of(primary), contigSupplementary());

        assertEquals(1, emitted.size());
        assertEquals("51M100S", emitted.get(0).getCigarString());
    }

    // The snap tests above run with no reference, which leaves the extension pass inert. These two supply one, so the
    // snap and the extension are exercised against each other.

    @Test
    public void testBoundarySnapSurvivesSoftClipExtension()
    {
        // The base a snap retracts is one bwa had aligned, so it still matches the reference and an extension running
        // afterwards reclaims it, putting the boundary back inside the intron. The read matches for its first 51 bases
        // and mismatches beyond, so that single base is the only one the extension could take.
        RefGenomeInterface refGenome = new TestGenome().with(CHR_1, 700, 'A').asRefGenome();
        byte[] readBases = bases("A".repeat(51) + "T".repeat(100));

        SAMRecord primary = primaryRecord(CHR_1, 150, "51M100S");
        primary.setReadBases(readBases);
        SAMRecord supp = supplementaryRecord(CHR_1, 499, "100S51M", CHR_1 + ",150,+,51M100S,60,0;");
        supp.setReadBases(readBases);

        List<SAMRecord> emitted = process(List.of(primary, supp), contigSupplementary(), null, refGenome, null);

        assertEquals(2, emitted.size());
        assertEquals("retraction not re-consumed by the extension", "50M101S", emitted.get(0).getCigarString());
        assertEquals(150, emitted.get(0).getAlignmentStart());
    }

    @Test
    public void testSoftClipExtensionStandsDownTheAlignmentScoreFloor()
    {
        // AS 20 is under the floor, but the extension lengthens 51M to 61M, so bwa's recorded score is stale and
        // pessimistic for the same reason a merge makes it stale. Split-read shape whose read coverage leaves a gap, so
        // no merge happens and the extension is the only thing that can stand the floor down.
        RefGenomeInterface refGenome = new TestGenome().with(CHR_1, 700, 'A').asRefGenome();
        byte[] readBases = bases("A".repeat(61) + "T".repeat(90));

        SAMRecord primary = primaryRecord(CHR_1, 150, "51M100S");
        primary.setReadBases(readBases);
        primary.setAttribute("AS", 20);
        SAMRecord supp = supplementaryRecord(CHR_1, 499, "100S51M", CHR_1 + ",150,+,51M100S,60,0;");
        supp.setReadBases(readBases);

        List<SAMRecord> emitted = process(List.of(primary, supp), contigSupplementary(), null, refGenome, null);

        assertEquals(2, emitted.size());
        assertFalse("an extended primary is not unmapped on a stale AS", emitted.get(0).getReadUnmappedFlag());
        assertEquals("61M90S", emitted.get(0).getCigarString());
    }
}
