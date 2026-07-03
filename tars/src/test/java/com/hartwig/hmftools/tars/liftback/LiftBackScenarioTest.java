package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.ReadRole.PRIMARY;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.mate;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.primary;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.supp;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

// Demonstrates the parametric LiftBackScenario harness: each test declares reads + tags + geometry and asserts the
// full pass-through result. Outcomes mirror the hand-wired LiftBackEndToEndTest so they are known-good.
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
        scenario()
                .read(primary("frag1", TX_CONTIG, 51, "100M"))
                .read(mate("frag1", TX_CONTIG, 197, "50M"))
                .run()
                .assertLifted("frag1", PRIMARY, CHR_1, 150, "50M100N50M");
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
                .assertSuppCount("frag3", 0);
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
