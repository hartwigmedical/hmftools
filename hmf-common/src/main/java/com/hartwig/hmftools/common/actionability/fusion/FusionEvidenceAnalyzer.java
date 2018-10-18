package com.hartwig.hmftools.common.actionability.fusion;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionEvidenceAnalyzer {

    @NotNull
    private final List<ActionableFusion> fusionPairs;
    @NotNull
    private final List<ActionablePromiscuousFive> promiscuousFive;
    @NotNull
    private final List<ActionablePromiscuousThree> promiscuousThree;

    FusionEvidenceAnalyzer(@NotNull final List<ActionableFusion> fusionPairs,
            @NotNull final List<ActionablePromiscuousFive> promiscuousFive,
            @NotNull final List<ActionablePromiscuousThree> promiscuousThree) {
        this.fusionPairs = fusionPairs;
        this.promiscuousFive = promiscuousFive;
        this.promiscuousThree = promiscuousThree;
    }

    @NotNull
    public Set<String> actionableGenes() {
        Set<String> genes = Sets.newHashSet();
        for (ActionableFusion fusionPairsSet : fusionPairs) {
            genes.add(fusionPairsSet.fiveGene());
        }
        for (ActionablePromiscuousThree promiscuousThreeSet : promiscuousThree) {
            genes.add(promiscuousThreeSet.gene());
        }
        for (ActionablePromiscuousFive promiscuousFiveSet : promiscuousFive) {
            genes.add(promiscuousFiveSet.gene());
        }
        return genes;
    }

    @NotNull
    public List<EvidenceItem> actionableFusions(@Nullable String doidsPrimaryTumorLocation,
            @NotNull CancerTypeAnalyzer cancerTypeAnalyzer, @NotNull GeneFusion geneFusion) {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        return evidenceItems;
    }
}
