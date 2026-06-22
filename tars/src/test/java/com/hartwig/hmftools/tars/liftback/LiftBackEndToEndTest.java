package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.assertLifted;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.secondMateRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.tars.common.TarsConstants;
import com.hartwig.hmftools.tars.liftback.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.RescueConfig;
import com.hartwig.hmftools.tars.liftback.tailextend.SoftclipTailExtender;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionConfig;
import com.hartwig.hmftools.tars.liftback.tailextend.TerminalMicroJunctionCollapser;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// End-to-end tests of how TarsApplication reacts to a bwa-mem2 read aligned against a transcript contig.
// Runs the FULL per-group pipeline: resolve -> junction-rescue -> terminal-collapse -> tail-extend ->
// canonicalize -> mate-patch -> NH/unmap policy. (Per-component behaviour is unit-tested in
// LiftBackResolverTest / SpliceLiftBackApplyTest / the rescue + tailextend suites.)
//
// Shares the standard three-exon contig, record builders and TestGenome with TarsTestFixtures (exons
// 100-199, 300-399, 500-549; introns 200-299, 400-499; contig length 250). REF-FREE passes (translation,
// SA/mate rewrite, NH, unmap-on-lift-failure) work with the default genome; REF-DEPENDENT passes (rescue
// ref-verify, tail-extend, collapse, canonicalize) only fire when the genome bases match the read / carry a
// splice motif -- set those explicitly via TestGenome.set() in the scenario.
public class LiftBackEndToEndTest
{
    // genomic introns implied by the exon layout, registered as annotated so rescue / tail-extend recognise them.
    private static Set<ChrBaseRegion> annotatedIntrons()
    {
        return Set.of(new ChrBaseRegion(CHR_1, 200, 299), new ChrBaseRegion(CHR_1, 400, 499));
    }

    // chr1 of all 'A' with canonical GT..AG seeded at the two intron boundaries, so the canonicalization /
    // rescue ref-verify passes have a motif to land on. Scenarios that need matching exon bases add them via set().
    private static TestGenome scenarioGenome(final int chr1Length)
    {
        return new TestGenome().with(CHR_1, chr1Length, 'A')
                .set(CHR_1, 200, "GT").set(CHR_1, 298, "AG")
                .set(CHR_1, 400, "GT").set(CHR_1, 498, "AG");
    }

    // set read bases (needed by the ref-dependent passes: tail-extend / collapse / canonicalize compare
    // the read's bases against the genome). Length should match the CIGAR's read-consuming length.
    private static SAMRecord withBases(final SAMRecord record, final String bases)
    {
        record.setReadBases(bases.getBytes());
        return record;
    }

    // ---- full-pipeline runner (mirrors LiftBackWorker's engine wiring) ----
    private static List<SAMRecord> runLiftBack(final TestGenome genome, final List<SAMRecord> reads)
    {
        RefSequenceSource ref = genome.asRefSource();
        AnnotatedJunctionIndex junctionIndex = new AnnotatedJunctionIndex(annotatedIntrons());

        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        JunctionRescueResolver rescue =
                new JunctionRescueResolver(junctionIndex, ref, RescueConfig.enabledDefaults());
        SoftclipTailExtender extender =
                new SoftclipTailExtender(ref, junctionIndex, TailExtensionConfig.enabledDefaults());
        TerminalMicroJunctionCollapser collapser =
                new TerminalMicroJunctionCollapser(ref, TarsConstants.MIN_JUNCTION_ANCHOR);
        JunctionCanonicalizer canonicalizer =
                new JunctionCanonicalizer(ref, TarsConstants.DEFAULT_MAX_SHIFT);

        LiftBackStats stats = new LiftBackStats();
        LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, rescue, extender, collapser, canonicalizer, ref, null, stats);

        List<SAMRecord> emitted = new ArrayList<>();
        processor.processNameGroup(reads, new LiftedMateInfoCache(), (record, result) -> emitted.add(record));
        return emitted;
    }

    // ============================ scenarios ============================

    @Test
    public void exonSpanningReadLiftsToJunctionCigar()
    {
        // R2 starts at tx 197 so its leading exon2 anchor is 4M (above the trimMicroAnchors threshold of 3).
        SAMRecord r1 = primaryRecord("read1", TX_CONTIG, 51, "100M");
        SAMRecord r2 = secondMateRecord("read1", TX_CONTIG, 197, "50M");

        runLiftBack(scenarioGenome(2000), List.of(r1, r2));

        // r1 spans exon1->exon2: 50M of exon1 + 100N intron1 + 50M exon2; mate fields point at the lifted R2.
        assertLifted(r1, CHR_1, 150, "50M100N50M");
        assertEquals(CHR_1, r1.getMateReferenceName());
        assertEquals(396, r1.getMateAlignmentStart());

        assertEquals(CHR_1, r2.getReferenceName());
        assertEquals(396, r2.getAlignmentStart());
        assertTrue(r2.getCigarString().contains("N"));
        assertEquals(150, r2.getMateAlignmentStart());
    }

    @Test
    public void unliftableReadIsMarkedUnmapped()
    {
        // starts past the contig end (tx 251 with contig len 250) -> translation fails -> emitted unmapped.
        SAMRecord r1 = primaryRecord("read2", TX_CONTIG, 251, "50M");
        SAMRecord r2 = secondMateRecord("read2", CHR_1, 600, "50M");

        runLiftBack(scenarioGenome(2000), List.of(r1, r2));

        assertTrue(r1.getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, r1.getCigarString());
        assertNull(r1.getStringAttribute("XA"));

        // the surviving mate loses proper-pair and is flagged mate-unmapped.
        assertEquals(CHR_1, r2.getReferenceName());
        assertFalse(r2.getProperPairFlag());
        assertTrue(r2.getMateUnmappedFlag());
    }

    @Test
    public void supplementarySaTagRewrittenToGenomicCoords()
    {
        // a split read: primary on exon1, supplementary on exon2, each pointing at the other via SA on the tx contig.
        SAMRecord r1 = primaryRecord("read4", TX_CONTIG, 51, "100M");
        r1.setAttribute("SA", TX_CONTIG + ",197,+,50M,60,0;");
        SAMRecord r1Supp = supplementaryRecord("read4", TX_CONTIG, 197, "50M");
        r1Supp.setAttribute("SA", TX_CONTIG + ",51,+,100M,60,0;");
        SAMRecord r2 = secondMateRecord("read4", CHR_1, 600, "50M");

        runLiftBack(scenarioGenome(2000), List.of(r1, r1Supp, r2));

        // SA must be lifted to chr1 and carry no _tx contig name.
        String r1Sa = r1.getStringAttribute("SA");
        assertTrue("SA should be lifted to chr1: " + r1Sa, r1Sa != null && r1Sa.startsWith(CHR_1 + ","));
        assertTrue("SA should not reference a _tx contig: " + r1Sa, !r1Sa.contains(TX_CONTIG));

        // supplementary lifts to its exon2 genomic locus, and both records' mate fields point at R2.
        assertEquals(CHR_1, r1Supp.getReferenceName());
        assertEquals(396, r1Supp.getAlignmentStart());
        assertEquals(CHR_1, r1.getMateReferenceName());
        assertEquals(600, r1.getMateAlignmentStart());
        assertEquals(CHR_1, r1Supp.getMateReferenceName());
        assertEquals(600, r1Supp.getMateAlignmentStart());
    }

    @Test
    public void terminalSoftclipExtendedIntoContiguousGenome()
    {
        // A read with a trailing soft-clip whose clipped bases continue contiguously in the genome (intron
        // retention, no junction). The tail-extender should walk the 10S into the reference -> 60M, with the
        // "tail-extended" note set. Genome is all 'A' away from the intron motifs, so all-'A' read bases match.
        TestGenome genome = scenarioGenome(2000);
        SAMRecord r1 = withBases(primaryRecord("read6", CHR_1, 1000, "50M10S"), "A".repeat(60));
        SAMRecord r2 = withBases(secondMateRecord("read6", CHR_1, 1500, "50M"), "A".repeat(50));

        runLiftBack(genome, List.of(r1, r2));

        assertEquals("60M", r1.getCigarString());
    }

    @Test
    public void splitReadRescuedAcrossAnnotatedJunction()
    {
        // bwa split the read into a primary (50M50S, exon1 side) and a supplementary (50S50M, exon2 side)
        // flanking the annotated intron 200-299. Rescue should merge them into one 50M100N50M primary and
        // drop the supp. The default genome already carries the canonical GT..AG at the intron boundaries,
        // so ref-verify passes.
        TestGenome genome = scenarioGenome(2000);
        SAMRecord primary = withBases(primaryRecord("read7", CHR_1, 150, "50M50S"), "A".repeat(100));
        SAMRecord supp = withBases(supplementaryRecord("read7", CHR_1, 300, "50S50M"), "A".repeat(100));
        SAMRecord mate = withBases(secondMateRecord("read7", CHR_1, 800, "50M"), "A".repeat(50));

        List<SAMRecord> emitted = runLiftBack(genome, List.of(primary, supp, mate));

        assertEquals("50M100N50M", primary.getCigarString());
        // the merged supplementary is dropped: only primary/1 + mate/2 remain.
        assertEquals(2, emitted.size());
    }

    @Test
    public void offMotifJunctionSlidToCanonical()
    {
        // A spliced read whose junction sits 2bp off a GT-AG motif. Canonicalization slides the intron +2 to
        // land on the motif: 10M10N10M -> 12M10N8M. Anchors are >= MIN_JUNCTION_ANCHOR so the terminal
        // collapser leaves the junction for canonicalization to handle (5M anchors would be collapsed first).
        TestGenome genome = scenarioGenome(2000)
                .set(CHR_1, 1, "CCCCCCCCCC")  // left exon, matches read[0..9]
                .set(CHR_1, 11, "CC")         // donor at shift 0 (non-canonical) + the +2 moved bases
                .set(CHR_1, 13, "GT")         // donor at shift +2 (canonical)
                .set(CHR_1, 21, "AG");        // acceptor at shift +2 (canonical)
        SAMRecord r1 = withBases(primaryRecord("read8", CHR_1, 1, "10M10N10M"), "C".repeat(20));
        SAMRecord r2 = withBases(secondMateRecord("read8", CHR_1, 500, "50M"), "A".repeat(50));

        runLiftBack(genome, List.of(r1, r2));

        assertEquals("12M10N8M", r1.getCigarString());
    }

    @Test
    public void nativeGenomicReadPassesThroughUnchanged()
    {
        SAMRecord r1 = primaryRecord("read3", CHR_1, 1000, "100M");
        SAMRecord r2 = secondMateRecord("read3", CHR_1, 1100, "50M");

        runLiftBack(scenarioGenome(2000), List.of(r1, r2));

        assertLifted(r1, CHR_1, 1000, "100M");
    }
}
