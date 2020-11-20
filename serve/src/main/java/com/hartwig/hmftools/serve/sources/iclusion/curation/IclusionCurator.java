package com.hartwig.hmftools.serve.sources.iclusion.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionCurator {

    private static final Logger LOGGER = LogManager.getLogger(IclusionCurator.class);

    @NotNull
    private final Set<String> evaluatedGenes = Sets.newHashSet();
    @NotNull
    private final Set<CurationEntry> evaluatedCurationEntries = Sets.newHashSet();

    public IclusionCurator() {
    }

    @NotNull
    public List<IclusionTrial> run(@NotNull List<IclusionTrial> trials) {
        List<IclusionTrial> curatedTrials = Lists.newArrayList();

        for (IclusionTrial trial : trials) {
            List<IclusionMutationCondition> curatedConditions = Lists.newArrayList();
            for (IclusionMutationCondition condition : trial.mutationConditions()) {
                List<IclusionMutation> curatedMutations = Lists.newArrayList();
                for (IclusionMutation mutation : condition.mutations()) {
                    curatedMutations.add(curate(mutation));
                }
                curatedConditions.add(ImmutableIclusionMutationCondition.builder().from(condition).mutations(curatedMutations).build());
            }
            curatedTrials.add(ImmutableIclusionTrial.builder().from(trial).mutationConditions(curatedConditions).build());
        }

        return curatedTrials;
    }

    public void reportUnusedCurationEntries() {
        int unusedEntryCount = 0;
        for (CurationEntry entry : CurationFactory.MUTATION_MAPPINGS.keySet()) {
            if (!evaluatedCurationEntries.contains(entry)) {
                unusedEntryCount++;
                LOGGER.warn("Entry '{}' hasn't been used during iClusion curation", entry);
            }
        }

        LOGGER.debug("Found {} unused iClusion curation entries. {} keys have been requested against {} entries",
                unusedEntryCount,
                evaluatedCurationEntries.size(),
                CurationFactory.MUTATION_MAPPINGS.size());

        int unusedGeneCount = 0;
        for (String gene : CurationFactory.GENE_MAPPINGS.keySet()) {
            if (!evaluatedGenes.contains(gene)) {
                unusedEntryCount++;
                LOGGER.warn("Gene '{}' hasn't been used during iClusion curation", gene);
            }
        }

        LOGGER.debug("Found {} unused iClusion curation genes. {} genes have been requested against {} genes",
                unusedGeneCount,
                evaluatedGenes.size(),
                CurationFactory.GENE_MAPPINGS.size());
    }

    @NotNull
    private IclusionMutation curate(@NotNull IclusionMutation mutation) {
        CurationEntry entry = new CurationEntry(mutation.gene(), mutation.name());
        evaluatedCurationEntries.add(entry);

        String gene = mutation.gene();
        evaluatedGenes.add(gene);

        IclusionMutation curatedMutation = mutation;
        if (CurationFactory.MUTATION_MAPPINGS.containsKey(entry)) {
            String mappedName = CurationFactory.MUTATION_MAPPINGS.get(entry).name();
            String mappedGene = CurationFactory.MUTATION_MAPPINGS.get(entry).gene();

            LOGGER.debug("Mapping mutation '{} on {}' to '{} on {}'", entry.name(), entry.gene(), mappedName, mappedGene);
            curatedMutation = ImmutableIclusionMutation.builder().from(curatedMutation).gene(mappedGene).name(mappedName).build();
        }

        if (CurationFactory.GENE_MAPPINGS.containsKey(gene)) {
            String mappedGene = CurationFactory.GENE_MAPPINGS.get(gene);
            LOGGER.debug("Mapping mutation '{} on {}' to '{} on {}'", entry.name(), entry.gene(), entry.name(), mappedGene);
            curatedMutation = ImmutableIclusionMutation.builder().from(curatedMutation).gene(mappedGene).build();
        }

        return curatedMutation;
    }
}
