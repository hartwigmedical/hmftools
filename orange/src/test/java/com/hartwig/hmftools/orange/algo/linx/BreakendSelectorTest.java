package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BreakendSelectorTest {

    @Test
    public void canSelectInterestingUnreportedBreakends() {
        LinxBreakend fiveGeneExon10 = createFive("five", "canonical", 10);
        LinxBreakend fiveGeneExon11 = createFive("five", "canonical", 11);
        LinxBreakend fiveGeneExon10NonCanonical = createFive("five", "non-canonical", 10);
        LinxBreakend fiveGeneExon20 = createFive("five", "canonical", 20);
        LinxBreakend threeGeneExon10 = createThree("three", "canonical", 10);
        LinxBreakend threeGeneExon20 = createThree("three", "canonical", 20);
        LinxBreakend threeGeneExon21 = createThree("three", "canonical", 21);
        LinxBreakend otherGene = createThree("other", "canonical", 10);
        List<LinxBreakend> allBreakends = Lists.newArrayList(fiveGeneExon10,
                fiveGeneExon11,
                fiveGeneExon10NonCanonical,
                fiveGeneExon20,
                threeGeneExon10,
                threeGeneExon20,
                threeGeneExon21,
                otherGene);

        KnownFusionData fiveData = knownFusion(KnownFusionType.PROMISCUOUS_5, "five", Strings.EMPTY);
        fiveData.setKnownExonData("canonical", "8;12", Strings.EMPTY);

        KnownFusionData threeData = knownFusion(KnownFusionType.PROMISCUOUS_3, Strings.EMPTY, "three");
        threeData.setKnownExonData("canonical", Strings.EMPTY, "20;21");

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.addData(fiveData);
        knownFusionCache.addData(threeData);

        LinxFusion fiveFusion = LinxTestFactory.fusionBuilder().geneStart("five").fusedExonUp(11).build();
        LinxFusion threeFusion = LinxTestFactory.fusionBuilder().geneEnd("three").fusedExonDown(21).build();
        List<LinxBreakend> potentiallyInteresting = BreakendSelector.selectInterestingUnreportedBreakends(allBreakends,
                Lists.newArrayList(fiveFusion, threeFusion),
                knownFusionCache);

        assertEquals(2, potentiallyInteresting.size());
        assertTrue(potentiallyInteresting.contains(fiveGeneExon10));
        assertTrue(potentiallyInteresting.contains(threeGeneExon20));
    }

    @NotNull
    private static LinxBreakend createFive(@NotNull String gene, @NotNull String transcript, int exon) {
        return create(gene, transcript, BreakendSelector.UPSTREAM_ORIENTATION, exon);
    }

    @NotNull
    private static LinxBreakend createThree(@NotNull String gene, @NotNull String transcript, int exon) {
        return create(gene, transcript, BreakendSelector.DOWNSTREAM_ORIENTATION, exon);
    }

    @NotNull
    private static LinxBreakend create(@NotNull String gene, @NotNull String transcript, @NotNull String geneOrientation, int exon) {
        return LinxTestFactory.breakendBuilder()
                .gene(gene)
                .transcriptId(transcript)
                .geneOrientation(geneOrientation)
                .nextSpliceExonRank(exon)
                .reportedDisruption(false)
                .disruptive(true)
                .build();
    }

    @NotNull
    private static KnownFusionData knownFusion(@NotNull KnownFusionType type, @NotNull String fiveGene, @NotNull String threeGene) {
        return new KnownFusionData(type, fiveGene, threeGene, Strings.EMPTY, Strings.EMPTY);
    }

}