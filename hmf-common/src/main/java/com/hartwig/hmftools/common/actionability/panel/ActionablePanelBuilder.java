package com.hartwig.hmftools.common.actionability.panel;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberType;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzer;

import org.jetbrains.annotations.NotNull;

public class ActionablePanelBuilder {

    private final Map<String, ImmutableActionablePanel.Builder> result = Maps.newHashMap();

    @NotNull
    public ActionablePanelBuilder addCopyNumbers(@NotNull final CopyNumberEvidenceAnalyzer copyNumberEvidenceAnalyzer) {
        copyNumberEvidenceAnalyzer.actionableGenes(CopyNumberType.DELETION).forEach(x -> select(x).deletion(true));
        copyNumberEvidenceAnalyzer.actionableGenes(CopyNumberType.AMPLIFICATION).forEach(x -> select(x).amplification(true));

        return this;
    }

    @NotNull
    public ActionablePanelBuilder addFusions(@NotNull final FusionEvidenceAnalyzer fusionAnalyser) {
        fusionAnalyser.actionableGenes().forEach(x -> select(x).fusion(true));
        return this;
    }

    @NotNull
    public ActionablePanelBuilder addVariants(@NotNull final SomaticVariantEvidenceAnalyzer variantEvidenceAnalyzer) {
        variantEvidenceAnalyzer.actionableGenes().forEach(x -> select(x).variant(true));
        return this;
    }

    @NotNull
    public List<ActionablePanel> build() {
        return result.values().stream().map(ImmutableActionablePanel.Builder::build).collect(Collectors.toList());
    }

    @NotNull
    private ImmutableActionablePanel.Builder select(@NotNull final String gene) {
        return result.computeIfAbsent(gene, this::create);
    }

    @NotNull
    private ImmutableActionablePanel.Builder create(@NotNull final String gene) {
        return ImmutableActionablePanel.builder().gene(gene).amplification(false).deletion(false).fusion(false).variant(false);
    }

}
