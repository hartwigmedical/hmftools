package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.SpliceCommon;
import com.hartwig.hmftools.tars.liftback.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.tars.liftback.rescue.ChrIntron;
import com.hartwig.hmftools.tars.liftback.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.RescueConfig;
import com.hartwig.hmftools.tars.liftback.tailextend.SoftclipTailExtender;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionConfig;
import com.hartwig.hmftools.tars.liftback.tailextend.TerminalMicroJunctionCollapser;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// End-to-end tests of how SpliceLiftBack reacts to a bwa-mem2 read aligned against a transcript contig.
// Unlike SpliceLiftBackEndToEndTest (which drives the resolver alone), this runs the FULL per-group
// pipeline: resolve -> junction-rescue -> terminal-collapse -> tail-extend -> canonicalize -> mate-patch
// -> NH/unmap policy, exactly as LiftBackWorker does in REDUX.
//
// REF-FREE passes (translation to genomic coords, SA/mate rewrite, NH, unmap-on-lift-failure) work with
// the default genome. REF-DEPENDENT passes (rescue ref-verify, tail-extend, collapse, canonicalize) only
// fire when the genome bases at the lifted coordinates match the read bases / carry a splice motif -- set
// those explicitly via Genome.set() in the scenario.
public class LiftBackScenarioTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    // Transcript contig (1-250) -> chr1. Exons: 100-199, 300-399, 500-549. Introns: 200-299, 400-499.
    private static final List<BaseRegion> EXONS =
            List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549));

    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1, EXONS);
    }

    // genomic introns implied by the exon layout, registered as annotated so rescue / tail-extend recognise them.
    private static Set<ChrIntron> annotatedIntrons()
    {
        return Set.of(new ChrIntron(CHR_1, 200, 299), new ChrIntron(CHR_1, 400, 499));
    }

    // ---- in-memory genome the scenario controls; coordinates are 1-based inclusive ----
    static final class Genome
    {
        private final Map<String,byte[]> mBases = new HashMap<>();

        Genome(final int chr1Length)
        {
            final byte[] seq = new byte[chr1Length];
            for(int i = 0; i < chr1Length; ++i)
                seq[i] = 'A';
            // place canonical donor (GT) / acceptor (AG) at the two intron boundaries so the
            // canonicalization pass has a motif to land on when a scenario exercises it.
            writeMotif(seq, 200, "GT");   // intron1 donor
            writeMotif(seq, 298, "AG");   // intron1 acceptor
            writeMotif(seq, 400, "GT");   // intron2 donor
            writeMotif(seq, 498, "AG");   // intron2 acceptor
            mBases.put(CHR_1, seq);
        }

        private static void writeMotif(final byte[] seq, final int oneBasedStart, final String motif)
        {
            for(int i = 0; i < motif.length(); ++i)
                seq[oneBasedStart - 1 + i] = (byte) motif.charAt(i);
        }

        // overwrite bases at a 1-based start with the given sequence (e.g. to match a read for NM=0 / extension)
        Genome set(final String chromosome, final int oneBasedStart, final String sequence)
        {
            final byte[] seq = mBases.get(chromosome);
            for(int i = 0; i < sequence.length(); ++i)
                seq[oneBasedStart - 1 + i] = (byte) sequence.charAt(i);
            return this;
        }

        RefSequenceSource asRefSource()
        {
            return TarsTestFixtures.refSource(CHR_1, new String(mBases.get(CHR_1), StandardCharsets.US_ASCII));
        }
    }

    // ---- read builders ----
    private static SAMRecord primary(
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

    // set read bases (needed by the ref-dependent passes: tail-extend / collapse / canonicalize compare
    // the read's bases against the genome). Length should match the CIGAR's read-consuming length.
    private static SAMRecord withBases(final SAMRecord record, final String bases)
    {
        record.setReadBases(bases.getBytes());
        return record;
    }

    private static SAMRecord supplementary(
            final String readName, final boolean firstOfPair, final String contig, final int pos, final String cigar)
    {
        final SAMRecord record = primary(readName, firstOfPair, contig, pos, cigar);
        record.setSupplementaryAlignmentFlag(true);
        return record;
    }

    // ---- full-pipeline runner (mirrors LiftBackWorker's engine wiring) ----
    private static List<SAMRecord> runLiftBack(final Genome genome, final List<SAMRecord> reads)
    {
        final RefSequenceSource ref = genome.asRefSource();
        final AnnotatedJunctionIndex junctionIndex = new AnnotatedJunctionIndex(annotatedIntrons());

        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        final JunctionRescueResolver rescue =
                new JunctionRescueResolver(junctionIndex, ref, RescueConfig.enabledDefaults());
        final SoftclipTailExtender extender =
                new SoftclipTailExtender(ref, junctionIndex, TailExtensionConfig.enabledDefaults());
        final TerminalMicroJunctionCollapser collapser =
                new TerminalMicroJunctionCollapser(ref, SpliceCommon.MIN_JUNCTION_ANCHOR);
        final JunctionCanonicalizer canonicalizer =
                new JunctionCanonicalizer(ref, JunctionCanonicalizer.DEFAULT_MAX_SHIFT);

        final LiftBackStats stats = new LiftBackStats();
        final LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, rescue, extender, collapser, canonicalizer, ref, 0, 0, stats);

        final List<SAMRecord> emitted = new ArrayList<>();
        processor.processNameGroup(reads, new LiftedMateInfoCache(), (record, result) -> emitted.add(record));
        return emitted;
    }

    // ============================ scenarios ============================

    @Test
    public void exonSpanningReadLiftsToJunctionCigar()
    {
        final SAMRecord r1 = primary("read1", true, TX_CONTIG, 51, "100M");
        final SAMRecord r2 = primary("read1", false, TX_CONTIG, 197, "50M");

        runLiftBack(new Genome(2000), List.of(r1, r2));

        // r1 spans exon1->exon2: 50M of exon1 + 100N intron1 + 50M exon2.
        assertEquals(CHR_1, r1.getReferenceName());
        assertEquals(150, r1.getAlignmentStart());
        assertEquals("50M100N50M", r1.getCigarString());
        assertTrue(r2.getCigarString().contains("N"));
    }

    @Test
    public void unliftableReadIsMarkedUnmapped()
    {
        // starts past the contig end (tx 251 with contig len 250) -> translation fails -> emitted unmapped.
        final SAMRecord r1 = primary("read2", true, TX_CONTIG, 251, "50M");
        final SAMRecord r2 = primary("read2", false, CHR_1, 600, "50M");

        runLiftBack(new Genome(2000), List.of(r1, r2));

        assertTrue(r1.getReadUnmappedFlag());
        assertTrue(r2.getMateUnmappedFlag());
    }

    @Test
    public void supplementarySaTagRewrittenToGenomicCoords()
    {
        // a split read: primary on exon1, supplementary on exon2, each pointing at the other via SA on the tx contig.
        final SAMRecord r1 = primary("read4", true, TX_CONTIG, 51, "100M");
        r1.setAttribute("SA", TX_CONTIG + ",197,+,50M,60,0;");
        final SAMRecord r1Supp = supplementary("read4", true, TX_CONTIG, 197, "50M");
        r1Supp.setAttribute("SA", TX_CONTIG + ",51,+,100M,60,0;");
        final SAMRecord r2 = primary("read4", false, CHR_1, 600, "50M");

        runLiftBack(new Genome(2000), List.of(r1, r1Supp, r2));

        // SA must be lifted to chr1 and carry no _tx contig name.
        final String r1Sa = r1.getStringAttribute("SA");
        assertTrue("SA should be lifted to chr1: " + r1Sa, r1Sa != null && r1Sa.startsWith(CHR_1 + ","));
        assertTrue("SA should not reference a _tx contig: " + r1Sa, !r1Sa.contains(TX_CONTIG));
    }

    @Test
    public void terminalSoftclipExtendedIntoContiguousGenome()
    {
        // A read with a trailing soft-clip whose clipped bases continue contiguously in the genome (intron
        // retention, no junction). The tail-extender should walk the 10S into the reference -> 60M, with the
        // "tail-extended" note set. Genome is all 'A' away from the intron motifs, so all-'A' read bases match.
        final Genome genome = new Genome(2000);
        final SAMRecord r1 = withBases(primary("read6", true, CHR_1, 1000, "50M10S"), "A".repeat(60));
        final SAMRecord r2 = withBases(primary("read6", false, CHR_1, 1500, "50M"), "A".repeat(50));

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
        final Genome genome = new Genome(2000);
        final SAMRecord primary = withBases(primary("read7", true, CHR_1, 150, "50M50S"), "A".repeat(100));
        final SAMRecord supp = withBases(supplementary("read7", true, CHR_1, 300, "50S50M"), "A".repeat(100));
        final SAMRecord mate = withBases(primary("read7", false, CHR_1, 800, "50M"), "A".repeat(50));

        final List<SAMRecord> emitted = runLiftBack(genome, List.of(primary, supp, mate));

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
        final Genome genome = new Genome(2000)
                .set(CHR_1, 1, "CCCCCCCCCC")  // left exon, matches read[0..9]
                .set(CHR_1, 11, "CC")         // donor at shift 0 (non-canonical) + the +2 moved bases
                .set(CHR_1, 13, "GT")         // donor at shift +2 (canonical)
                .set(CHR_1, 21, "AG");        // acceptor at shift +2 (canonical)
        final SAMRecord r1 = withBases(primary("read8", true, CHR_1, 1, "10M10N10M"), "C".repeat(20));
        final SAMRecord r2 = withBases(primary("read8", false, CHR_1, 500, "50M"), "A".repeat(50));

        runLiftBack(genome, List.of(r1, r2));

        assertEquals("12M10N8M", r1.getCigarString());
    }

    @Test
    public void nativeGenomicReadPassesThroughUnchanged()
    {
        final SAMRecord r1 = primary("read3", true, CHR_1, 1000, "100M");
        final SAMRecord r2 = primary("read3", false, CHR_1, 1100, "50M");

        runLiftBack(new Genome(2000), List.of(r1, r2));

        assertEquals(1000, r1.getAlignmentStart());
        assertEquals("100M", r1.getCigarString());
    }
}
