package com.hartwig.hmftools.common.actionability.panel;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.Actionable;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberType;
import com.hartwig.hmftools.common.actionability.drup.DrupActionabilityModel;
import com.hartwig.hmftools.common.actionability.drup.DrupActionabilityModelFactory;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzer;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

// NOTE: PLEASE DO NOT DELETE THIS... CURRENTLY USED AS PART OF THE GENE PANEL PROCESS
// TODO: Move this somewhere it can be run automatically, like into the knowledgebase perhaps?
@SuppressWarnings("unused")
public class ActionablePanelBuilder {

    private static final Set<String> EXCLUDE = Sets.newHashSet("iclusion");

    private final Map<String, ImmutableActionablePanel.Builder> result = Maps.newHashMap();

    public ActionablePanelBuilder(@NotNull final String knowledgebaseDirector, @NotNull final String drupLocation) throws IOException {
        final ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirector);
        final DrupActionabilityModel drup = DrupActionabilityModelFactory.buildFromCsv(drupLocation);

        addCopyNumbers(actionabilityAnalyzer.cnvAnalyzer());
        addVariants(actionabilityAnalyzer.variantAnalyzer());
        addDrup(drup);
    }

    @NotNull
    public ActionablePanelBuilder addCopyNumbers(@NotNull final CopyNumberEvidenceAnalyzer copyNumberEvidenceAnalyzer) {
        copyNumberEvidenceAnalyzer.actionableCopyNumbers().stream().filter(this::filterSource).forEach(x -> {
            final ImmutableActionablePanel.Builder builder = addActionable(x.gene(), x);
            if (x.type() == CopyNumberType.AMPLIFICATION) {
                builder.amplification(true);
            } else {
                builder.deletion(true);
            }
        });
        return this;
    }

    @NotNull
    public ActionablePanelBuilder addVariants(@NotNull final SomaticVariantEvidenceAnalyzer variantEvidenceAnalyzer) {
        variantEvidenceAnalyzer.actionableRanges()
                .stream()
                .filter(this::filterSource)
                .forEach(x -> addActionable(x.gene(), x).variant(true));
        variantEvidenceAnalyzer.actionableVariants()
                .stream()
                .filter(this::filterSource)
                .forEach(x -> addActionable(x.gene(), x).variant(true));
        return this;
    }

    @NotNull
    public ActionablePanelBuilder addDrup(@NotNull final DrupActionabilityModel drup) {
        drup.actionableGenes().forEach(this::addDrup);
        return this;
    }

    @NotNull
    public List<ActionablePanel> build() {
        return result.values().stream().map(ImmutableActionablePanel.Builder::build).collect(Collectors.toList());
    }

    @NotNull
    private ImmutableActionablePanel.Builder addActionable(@NotNull final String gene, @NotNull final Actionable actionable) {
        final ImmutableActionablePanel.Builder builder = select(gene);
        final ActionablePanel current = builder.build();

        if (actionable.response().equals("Responsive")) {
            boolean currentEmpty = current.responsive().isEmpty();
            boolean higherLevel = actionable.level().compareTo(current.responsive()) < 0;
            if (currentEmpty || higherLevel) {
                builder.responsive(actionable.level()).responsiveSource(actionable.source());
            }
        } else {
            boolean currentEmpty = current.resistant().isEmpty();
            boolean higherLevel = actionable.level().compareTo(current.resistant()) < 0;
            if (currentEmpty || higherLevel) {
                builder.resistant(actionable.level()).resistantSource(actionable.source());
            }
        }

        return builder;
    }

    @NotNull
    private ActionablePanelBuilder addDrup(@NotNull final String gene) {
        select(gene).drup(true);
        return this;
    }

    @NotNull
    private ImmutableActionablePanel.Builder select(@NotNull final String gene) {
        return result.computeIfAbsent(gene, this::create);
    }

    @NotNull
    private ImmutableActionablePanel.Builder create(@NotNull final String gene) {
        return ImmutableActionablePanel.builder()
                .gene(gene)
                .amplification(false)
                .deletion(false)
                .variant(false)
                .drup(false)
                .responsive(Strings.EMPTY)
                .responsiveSource(Strings.EMPTY)
                .resistant(Strings.EMPTY)
                .resistantSource(Strings.EMPTY);
    }

    private boolean filterSource(@NotNull final Actionable actionable) {
        return !EXCLUDE.contains(actionable.source());
    }

}
