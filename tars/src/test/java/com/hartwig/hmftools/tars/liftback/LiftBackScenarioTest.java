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

// Full-pipeline scenarios from candidate construction through BAM emission.
public class LiftBackScenarioTest
{
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
    public void testExonSpanningReadLiftsToJunctionCigar()
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
    public void testUnliftableReadIsMarkedUnmapped()
    {
        // primary starts past the contig end (tx 251, contig length 250) -> translation fails -> emitted unmapped
        // beside its mapped mate; the surviving mate loses proper-pair and is flagged mate-unmapped.
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
    public void testSupplementarySaTagRewrittenToGenomicCoords()
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
    public void testGenomicTerminalSoftclipIsPreserved()
    {
        // TARS does not opportunistically turn a terminal soft clip into matches.
        scenario()
                .read(primary("frag9", CHR_1, 1000, "50M10S").bases("A".repeat(60)))
                .read(mate("frag9", CHR_1, 1500, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag9", PRIMARY, CHR_1, 1000, "50M10S");
    }

    @Test
    public void testTxTerminalSoftclipIsPreserved()
    {
        // Liftback preserves a terminal soft clip unless a documented overhang collapse or supplementary merge changes it.
        scenario()
                .read(primary("frag10", TX_CONTIG, 1, "40M10S").bases("A".repeat(50)))
                .read(mate("frag10", CHR_1, 1500, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag10", PRIMARY, CHR_1, 100, "40M10S");
    }

    @Test
    public void testSingleLocusTxMapQualityZeroBumpsToSixty()
    {
        // a tx-only read at a single genomic locus with bwa MAPQ 0 lifts and is promoted to 60 with no XA alt.
        scenario()
                .read(primary("frag2", TX_CONTIG, 1, "50M").mapQuality(0))
                .run()
                .assertLifted("frag2", PRIMARY, CHR_1, 100, "50M")
                .assertMapQuality("frag2", PRIMARY, 60)
                .assertNoXa("frag2", PRIMARY);
    }

    @Test
    public void testSplitReadResolvedAcrossAnnotatedJunctionDropsSupp()
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
    public void testSplitReadUniquePairBumpsMapQualityToSixty()
    {
        // bwa gave the split halves MAPQ 0; the XA alt overlaps the primary placement, so the primary+supplementary
        // pair maps to a single locus and the merged spliced primary is promoted to 60.
        scenario()
                .read(primary("frag5", CHR_1, 150, "50M50S").mapQuality(0).xa(CHR_1 + ",+150,50M50S,0").bases("A".repeat(100)))
                .read(supp("frag5", CHR_1, 300, "50S50M").mapQuality(0).bases("A".repeat(100)))
                .read(mate("frag5", CHR_1, 800, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag5", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertSuppCount("frag5", 0)
                .assertMapQuality("frag5", PRIMARY, 60);
    }

    @Test
    public void testSupplementaryCanMergeWithXaCandidateBeforePrimaryChoice()
    {
        // The primary's XA candidate and supplementary form the supported splice. The reciprocal SA tags match BWA output.
        scenario()
                .read(primary("frag7", CHR_1, 900, "50M50S").mapQuality(0)
                        .xa(CHR_1 + ",+150,50M50S,0")
                        .sa(CHR_1 + ",300,+,50S50M,0,0;")
                        .bases("A".repeat(100)))
                .read(supp("frag7", CHR_1, 300, "50S50M").mapQuality(0)
                        .sa(CHR_1 + ",900,+,50M50S,0,0;")
                        .bases("A".repeat(100)))
                .read(mate("frag7", CHR_1, 800, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag7", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertSuppCount("frag7", 0)
                .assertMapQuality("frag7", PRIMARY, 0)
                .assertXa("frag7", PRIMARY, CHR_1 + ",+900,50M50S,0;");
    }

    @Test
    public void testSplitReadMultiLocusPairKeepsBwaMapQuality()
    {
        // same merge, but the XA alt is a distinct locus (chr1:900), so the pair does not uniquely map: the merged
        // primary keeps its bwa MAPQ 0 rather than being bumped.
        scenario()
                .read(primary("frag6", CHR_1, 150, "50M50S").mapQuality(0).xa(CHR_1 + ",+900,50M50S,0").bases("A".repeat(100)))
                .read(supp("frag6", CHR_1, 300, "50S50M").mapQuality(0).bases("A".repeat(100)))
                .read(mate("frag6", CHR_1, 800, "50M").bases("A".repeat(50)))
                .run()
                .assertLifted("frag6", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertSuppCount("frag6", 0)
                .assertMapQuality("frag6", PRIMARY, 0);
    }

    @Test
    public void testNativeGenomicReadPassesThroughUntouched()
    {
        scenario()
                .read(primary("frag4", CHR_1, 1000, "100M"))
                .read(mate("frag4", CHR_1, 1100, "50M"))
                .run()
                .assertLifted("frag4", PRIMARY, CHR_1, 1000, "100M");
    }
}
