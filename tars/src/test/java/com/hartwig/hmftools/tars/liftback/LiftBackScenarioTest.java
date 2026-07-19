package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.ReadRole.MATE;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.ReadRole.PRIMARY;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.ReadRole.SUPPLEMENTARY;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.mate;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.primary;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.supp;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

// Whole-pipeline tests driven through the parametric LiftBackScenario harness: each test declares reads + tags +
// geometry and asserts the lifted result of the full per-group pipeline (resolve -> overhang-gate peel ->
// supplementary-resolve -> post-resolve peel + reclaim -> mate-patch -> NH/unmap policy). Per-component behaviour is
// unit-tested in LiftBackResolverTest / SpliceLiftBackApplyTest / the supplementary + overhang suites.
//
// Standard three-exon contig (exons 100-199, 300-399, 500-549; introns 200-299, 400-499; contig length 250) with an
// all-'A' genome carrying canonical GT..AG at the intron boundaries so the ref-dependent passes (supplementary-resolve
// ref-verify, tail-extend, collapse) have a motif to land on; scenarios that need matching exon bases add them via bases().
public class LiftBackScenarioTest
{
    // three-exon contig geometry (exons 100-199, 300-399, 500-549; introns 200-299, 400-499), all-'A' genome with
    // canonical GT..AG seeded at the intron boundaries so supplementary-resolve ref-verify / canonicalize have a motif to land on.
    private static TestGenome scenarioGenome()
    {
        return new TestGenome().with(CHR_1, 2000, 'A')
                .set(CHR_1, 200, "GT").set(CHR_1, 298, "AG")
                .set(CHR_1, 400, "GT").set(CHR_1, 498, "AG");
    }

    private static LiftBackScenario scenario()
    {
        return LiftBackScenario.create()
                .contig(threeExonContig())
                .genome(scenarioGenome())
                .annotatedIntron(CHR_1, 200, 299)
                .annotatedIntron(CHR_1, 400, 499);
    }

    @Test
    public void exonSpanningReadLiftsToJunctionCigar()
    {
        // r1 spans exon1->exon2 (50M100N50M); the mate at tx 197 crosses exon2->exon3 with a 4M leading anchor
        // (above the micro-anchor floor of 3) -> 4M100N46M. Mate fields cross-point at each other's lifted locus.
        scenario()
                .read(primary("frag1", TX_CONTIG, 51, "100M"))
                .read(mate("frag1", TX_CONTIG, 197, "50M"))
                .run()
                .assertLifted("frag1", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertLifted("frag1", MATE, CHR_1, 396, "4M100N46M")
                .assertMate("frag1", PRIMARY, CHR_1, 396)
                .assertMate("frag1", MATE, CHR_1, 150);
    }

    @Test
    public void unliftableReadIsMarkedUnmapped()
    {
        // primary starts past the contig end (tx 251, contig length 250) -> translation fails -> emitted unmapped,
        // XA stripped; the surviving mate loses proper-pair and is flagged mate-unmapped.
        scenario()
                .read(primary("frag7", TX_CONTIG, 251, "50M"))
                .read(mate("frag7", CHR_1, 600, "50M"))
                .run()
                .assertUnmapped("frag7", PRIMARY)
                .assertNoXa("frag7", PRIMARY)
                .assertLifted("frag7", MATE, CHR_1, 600, "50M")
                .assertProperPair("frag7", MATE, false)
                .assertMateUnmapped("frag7", MATE);
    }

    @Test
    public void supplementarySaTagRewrittenToGenomicCoords()
    {
        // a split read: primary on exon1 (100M, no terminal softclip so no merge), supplementary on exon2, each
        // pointing at the other via a tx-contig SA. The primary's SA must be lifted to chr1 (no _tx name), the supp
        // lifts to its exon2/exon3 genomic locus, and both records' mate fields point at the mate.
        scenario()
                .read(primary("frag8", TX_CONTIG, 51, "100M").sa(TX_CONTIG + ",197,+,50M,60,0;"))
                .read(supp("frag8", TX_CONTIG, 197, "50M").sa(TX_CONTIG + ",51,+,100M,60,0;"))
                .read(mate("frag8", CHR_1, 600, "50M"))
                .run()
                .assertSaGenomic("frag8", PRIMARY, CHR_1)
                .assertLifted("frag8", SUPPLEMENTARY, CHR_1, 396, "4M100N46M")
                .assertMate("frag8", PRIMARY, CHR_1, 600)
                .assertMate("frag8", SUPPLEMENTARY, CHR_1, 600);
    }

    @Test
    public void genomicTerminalSoftclipNotReclaimed()
    {
        // a GENOMIC (non-tx) read with a trailing soft-clip whose clipped bases continue contiguously in the all-'A'
        // genome. The standalone reclaim is tx-match-only, so the genomic over-clip is left as bwa placed it: 50M10S stays.
        scenario()
                .read(primary("frag9", CHR_1, 1000, "50M10S").bases("A".repeat(60)))
                .read(mate("frag9", CHR_1, 1500, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag9", PRIMARY, CHR_1, 1000, "50M10S");
    }

    @Test
    public void txMatchTerminalSoftclipReclaimed()
    {
        // a tx-contig read (exon1, lifts to chr1:100) with a trailing soft-clip whose clipped bases continue
        // contiguously in the genome. The chosen primary is a tx-match, so the post-resolve reclaim walks the 10S into
        // the reference -> 50M.
        scenario()
                .read(primary("frag10", TX_CONTIG, 1, "40M10S").bases("A".repeat(50)))
                .read(mate("frag10", CHR_1, 1500, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag10", PRIMARY, CHR_1, 100, "50M");
    }

    @Test
    public void singleLocusTxMapqZeroBumpsToSixty()
    {
        // a tx-only read at a single genomic locus with bwa MAPQ 0 lifts and is promoted to 60 with no XA alt.
        scenario()
                .read(primary("frag2", TX_CONTIG, 1, "50M").mapq(0))
                .run()
                .assertLifted("frag2", PRIMARY, CHR_1, 100, "50M")
                .assertMapq("frag2", PRIMARY, 60)
                .assertNoXa("frag2", PRIMARY);
    }

    @Test
    public void splitReadResolvedAcrossAnnotatedJunctionDropsSupp()
    {
        // primary (exon1 side) + supplementary (exon2 side) flanking the annotated intron 200-299 merge into one
        // spliced primary; the supplementary is dropped, leaving primary/1 + mate/2.
        scenario()
                .read(primary("frag3", CHR_1, 150, "50M50S").bases("A".repeat(100)))
                .read(supp("frag3", CHR_1, 300, "50S50M").bases("A".repeat(100)))
                .read(mate("frag3", CHR_1, 800, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag3", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertSuppCount("frag3", 0)
                .assertEmittedCount(2); // merged supp dropped: only primary/1 + mate/2 remain
    }

    @Test
    public void splitReadUniquePairBumpsMapqToSixty()
    {
        // bwa gave the split halves MAPQ 0; the XA alt overlaps the primary placement, so the primary+supplementary
        // pair maps to a single locus and the merged spliced primary is promoted to 60.
        scenario()
                .read(primary("frag5", CHR_1, 150, "50M50S").mapq(0).xa(CHR_1 + ",+150,50M50S,0").bases("A".repeat(100)))
                .read(supp("frag5", CHR_1, 300, "50S50M").mapq(0).bases("A".repeat(100)))
                .read(mate("frag5", CHR_1, 800, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag5", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertSuppCount("frag5", 0)
                .assertMapq("frag5", PRIMARY, 60);
    }

    @Test
    public void splitReadMultiLocusPairKeepsBwaMapq()
    {
        // same merge, but the XA alt is a distinct locus (chr1:900), so the pair does not uniquely map: the merged
        // primary keeps its bwa MAPQ 0 rather than being bumped.
        scenario()
                .read(primary("frag6", CHR_1, 150, "50M50S").mapq(0).xa(CHR_1 + ",+900,50M50S,0").bases("A".repeat(100)))
                .read(supp("frag6", CHR_1, 300, "50S50M").mapq(0).bases("A".repeat(100)))
                .read(mate("frag6", CHR_1, 800, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag6", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertSuppCount("frag6", 0)
                .assertMapq("frag6", PRIMARY, 0);
    }

    @Test
    public void nativeGenomicReadPassesThroughUntouched()
    {
        scenario()
                .read(primary("frag4", CHR_1, 1000, "100M"))
                .read(mate("frag4", CHR_1, 1100, "50M"))
                .run()
                .assertLifted("frag4", PRIMARY, CHR_1, 1000, "100M");
    }
}
