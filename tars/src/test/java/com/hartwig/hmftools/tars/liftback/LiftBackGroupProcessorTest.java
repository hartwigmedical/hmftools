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
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryConfig;
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
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// A record that stays supplementary must carry an SA tag, since REDUX dedup reads the primary's coords from it: a supp whose
// SA entries all fail to lift is dropped rather than emitted with a null SA.
public class LiftBackGroupProcessorTest
{
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

    // a non-null resolver activates the AS-unmap gate, but with no annotated junctions it cannot improve any primary
    private static SupplementaryResolver noopSupplementary()
    {
        return new SupplementaryResolver(Collections.emptySet(), supplementaryConfig());
    }

    @Test
    public void testPrimaryLiftingIntoExcludedRegionIsUnmapped()
    {
        // the tx primary lifts to chr1:100, inside the excluded region; REDUX-style unmapping keeps the record but clears the cigar
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
        // the genomic ref stub is all 'A' and the read ends in two mismatches, so NM must be recomputed to 2 rather than carried
        // from the stale tx-contig NM:0, and MD must be dropped
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
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        assertEquals(Integer.valueOf(1), emitted.get(0).getIntegerAttribute("NH"));
    }

    @Test
    public void testOverCapGenomicPrimaryMapQuality0NoXaUnmapped()
    {
        // bwa emits MAPQ 0 and no XA when a read maps past the -h 75 XA cap, so this genomic primary is unmapped even though the
        // discriminator sees a single locus and would otherwise bump it to 60
        SAMRecord primary = primaryRecord(CHR_1, 100, "50M");
        primary.setMappingQuality(0);

        List<SAMRecord> emitted = process(List.of(primary));

        assertEquals(1, emitted.size());
        assertTrue(emitted.get(0).getReadUnmappedFlag());
    }

    @Test
    public void testOverCapTxPrimaryMapQuality0NoXaKept()
    {
        // a tx-contig primary at MAPQ 0 with no XA hit the XA cap on transcript contigs of one gene which all lift to one genomic
        // locus, so the over-cap rule is ref-only-gated and must not unmap it
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
        // MAPQ 0 with XA present is an ordinary multimapper within the cap: keep and lift it
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
        // /1 lifts to chr1:100 and /2 (exon2) to chr1:300; the per-group pair must point each mate's fields at the other's lifted
        // coords without a whole-sample first pass
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
        // /1 is unmapped by the genomic over-cap rule while /2 stays mapped. REDUX parks the unmapped mate at /2's genomic
        // coordinate, and /2 points its unmapped-mate fields back at itself.
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
        // Ultima single-end: an unpaired primary lifts with no mate patching, and getFirstOfPairFlag() throws on unpaired records
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
        // exercises the mate refresh with a supp in the group, which is where the unpaired getFirstOfPairFlag() call sat
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
        // bwa can emit the same junction on multiple tx contigs, so supps lifting to identical coords and strand collapse to one
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
        // the supp's only SA entry is an out-of-range tx position, so the SA rewrite yields null and the supp is dropped
        SAMRecord primary = primaryRecord(TX_CONTIG, 100, "50M");
        SAMRecord supp = supplementaryRecord(TX_CONTIG, 110, "30M", TX_CONTIG + ",9999,+,30M,0,0;");

        List<SAMRecord> emitted = process(List.of(primary, supp));

        assertEquals(1, emitted.size());
        assertFalse(emitted.get(0).getSupplementaryAlignmentFlag());
    }

    @Test
    public void testSupplementaryWithLiftableSaIsKept()
    {
        // same shape as the orphan case, but the SA entry lifts, so the supp is emitted with a rewritten genomic SA
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
        // the primary lifts but scores below the AS floor of 30, and supplementary resolve with no junctions cannot improve it
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
        // the AS-floor unmap happens during decide, so the supp-orphan gate sees the unmapped primary and drops the supplementary.
        // Unmapping only at emit would leave the supp with an SA referencing a primary emitted as unmapped, a dangling locus for
        // REDUX FragmentCoords. The supp's SA lifts cleanly, so without the decide-time unmap this emits two records, not one.
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
        // a fully-overlapping short-insert pair whose mates map equally well (NM 0, tied genomic score) to two genomic copies:
        // the fragment must not be torn across the copies, so both mates land on one chromosome and the other goes to XA
        String seq = "ACGTTGCA".repeat(7).substring(0, 55);
        RefGenomeInterface ref = new TestGenome()
                .with("chr5", 300, 'A').with("chr10", 300, 'A')
                .set("chr5", 100, seq).set("chr10", 100, seq)
                .asRefGenome();

        SAMRecord mate1 = primaryRecord("frag", "chr5", 100, "55M");
        mate1.setMappingQuality(0);
        mate1.setReadNegativeStrandFlag(true);
        mate1.setProperPairFlag(false);
        mate1.setReadBases(bases(seq));
        mate1.setAttribute("XA", "chr10,-100,55M,0;");

        SAMRecord mate2 = secondMateRecord("frag", "chr10", 100, "55M");
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
                supplementaryConfig());

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

    // bwa emits a pair with both ends unmapped as flags 77/141. Clearing the mate-unmapped bit leaves a record claiming a mapped
    // mate with RNEXT unset, which REDUX rejects with htsjdk's INVALID_FLAG_MATE_UNMAPPED.
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

    // junctions built from the same sidecar entry the discriminator lifts against, so the intron coords and chromosome key match
    // what the lift emits: chr1 introns 200-299 and 400-499 between the three exons
    private static SupplementaryResolver contigSupplementary()
    {
        return new SupplementaryResolver(
                EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())), null, supplementaryConfig());
    }

    @Test
    public void testOverhangCollapseAlsoRunsOnASupplementary()
    {
        // a weak 2M overhang before a supplementary's splice junction must be collapsed like a primary's. Left ungated the 2M
        // block is the fragment's lowest mapped coordinate and becomes the reported fusion junction, two bases of evidence
        // standing in for the real anchor past the N gap.
        RefGenomeInterface ref = new TestGenome()
                .with(CHR_1, 6000, 'A')
                .asRefGenome();

        // the primary's matched bases stop at 50 and the supplementary's resume at 100, so the read-coverage gap refuses the merge
        // and the supplementary is emitted on its own, which is where an ungated overhang survives
        SAMRecord primary = primaryRecord(CHR_1, 5000, "50M101S");
        primary.setReadBases(bases("A".repeat(151)));

        SAMRecord supp = supplementaryRecord(CHR_1, 1000, "99S2M100N50M", CHR_1 + ",5000,+,50M101S,60,0;");
        supp.setReadBases(primary.getReadBases());

        List<SAMRecord> emitted = process(List.of(primary, supp), contigSupplementary(), new OverhangGate(ref), ref, null);

        SAMRecord emittedPrimary = emitted.stream().filter(x -> !x.getSupplementaryAlignmentFlag()).findFirst().orElse(null);
        SAMRecord emittedSupp = emitted.stream().filter(SAMRecord::getSupplementaryAlignmentFlag).findFirst().orElse(null);
        assertTrue(emittedPrimary != null);
        assertTrue(emittedSupp != null);
        assertEquals(5000, emittedPrimary.getAlignmentStart());
        assertEquals("50M101S", emittedPrimary.getCigarString());
        assertEquals(1001, emittedSupp.getAlignmentStart());
        assertEquals("151M", emittedSupp.getCigarString());
        assertEquals(Integer.valueOf(0), emittedPrimary.getIntegerAttribute("NM"));
        assertEquals(Integer.valueOf(0), emittedSupp.getIntegerAttribute("NM"));
        assertEquals(CHR_1 + ",1001,+,151M,60,0;", emittedPrimary.getStringAttribute("SA"));
        assertEquals(CHR_1 + ",5000,+,50M101S,60,0;", emittedSupp.getStringAttribute("SA"));
    }

    @Test
    public void testEveryEmittedAlignmentGetsACompleteReciprocalSaTag()
    {
        SAMRecord primary = primaryRecord(CHR_1, 100, "50M100S");
        SAMRecord supp1 = supplementaryRecord(
                CHR_1, 300, "50H50M50H", CHR_1 + ",100,+,50M100S,60,0;");
        SAMRecord supp2 = supplementaryRecord(
                CHR_1, 500, "100H50M", CHR_1 + ",100,+,50M100S,60,0;");
        supp2.setReadNegativeStrandFlag(true);
        supp2.setMappingQuality(40);

        List<SAMRecord> emitted = process(List.of(primary, supp1, supp2));

        assertEquals(3, emitted.size());
        SAMRecord emittedPrimary = emitted.get(0);
        SAMRecord emittedSupp1 = emitted.get(1);
        SAMRecord emittedSupp2 = emitted.get(2);
        assertEquals(
                CHR_1 + ",300,+,50S50M50S,60,0;" + CHR_1 + ",500,-,100S50M,40,0;",
                emittedPrimary.getStringAttribute("SA"));
        assertEquals(
                CHR_1 + ",100,+,50M100S,60,0;" + CHR_1 + ",500,-,100S50M,40,0;",
                emittedSupp1.getStringAttribute("SA"));
        assertEquals(
                CHR_1 + ",100,+,50M100S,60,0;" + CHR_1 + ",300,+,50S50M50S,60,0;",
                emittedSupp2.getStringAttribute("SA"));
    }

}
