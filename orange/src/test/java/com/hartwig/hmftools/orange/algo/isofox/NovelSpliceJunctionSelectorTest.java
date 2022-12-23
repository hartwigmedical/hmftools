package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCacheTestFactory;
import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class NovelSpliceJunctionSelectorTest {

    @Test
    public void canSelectSkippedExons() {
        NovelSpliceJunction match = create("gene 1", AltSpliceJunctionType.SKIPPED_EXONS, 10, 1);
        NovelSpliceJunction tooFewFragments = create("gene 1", AltSpliceJunctionType.SKIPPED_EXONS, 2, 1);
        NovelSpliceJunction tooHighCohortFreq = create("gene 1", AltSpliceJunctionType.SKIPPED_EXONS, 20, 50);
        NovelSpliceJunction alreadyHasFusion = create("gene 2", AltSpliceJunctionType.SKIPPED_EXONS, 20, 2);
        NovelSpliceJunction noKnownFusion = create("gene 3", AltSpliceJunctionType.SKIPPED_EXONS, 20, 2);
        List<NovelSpliceJunction> junctions =
                Lists.newArrayList(match, tooFewFragments, tooHighCohortFreq, alreadyHasFusion, noKnownFusion);

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.addData(KnownFusionCacheTestFactory.createExonDelDup("gene 1"));
        knownFusionCache.addData(KnownFusionCacheTestFactory.createExonDelDup("gene 2"));

        List<LinxFusion> linxFusions = Lists.newArrayList(LinxTestFactory.fusionBuilder().geneStart("gene 2").geneEnd("gene 2").build());

        List<NovelSpliceJunction> skippedExons = NovelSpliceJunctionSelector.selectSkippedExons(junctions, linxFusions, knownFusionCache);
        assertEquals(1, skippedExons.size());
        assertEquals(match, skippedExons.get(0));
    }

    @Test
    public void canSelectNovelExonIntrons() {
        NovelSpliceJunction match = create("gene 1", AltSpliceJunctionType.NOVEL_INTRON, 10, 1);
        NovelSpliceJunction tooFewFragments = create("gene 1", AltSpliceJunctionType.NOVEL_INTRON, 2, 1);
        NovelSpliceJunction tooHighCohortFreq = create("gene 1", AltSpliceJunctionType.NOVEL_INTRON, 10, 20);
        NovelSpliceJunction noDriverGene = create("gene 2", AltSpliceJunctionType.NOVEL_INTRON, 10, 1);
        List<NovelSpliceJunction> junctions = Lists.newArrayList(match, tooFewFragments, tooHighCohortFreq, noDriverGene);

        List<DriverGene> driverGenes = Lists.newArrayList(DriverGeneTestFactory.builder().gene("gene 1").build());

        List<NovelSpliceJunction> novelExonsIntrons = NovelSpliceJunctionSelector.selectNovelExonsIntrons(junctions, driverGenes);
        assertEquals(1, novelExonsIntrons.size());
        assertEquals(match, novelExonsIntrons.get(0));
    }

    @NotNull
    private static NovelSpliceJunction create(@NotNull String gene, @NotNull AltSpliceJunctionType type, int fragmentCount,
            int cohortFrequency) {
        return IsofoxTestFactory.novelSpliceJunctionBuilder()
                .geneName(gene)
                .type(type)
                .fragmentCount(fragmentCount)
                .cohortFrequency(cohortFrequency)
                .build();
    }

}