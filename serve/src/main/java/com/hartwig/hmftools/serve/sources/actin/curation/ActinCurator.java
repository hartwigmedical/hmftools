package com.hartwig.hmftools.serve.sources.actin.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ActinCurator {

    private static final Logger LOGGER = LogManager.getLogger(ActinCurator.class);

    @NotNull
    private final Set<String> evaluatedGenes = Sets.newHashSet();
    @NotNull
    private final Set<CurationEntry> evaluatedCurationEntries = Sets.newHashSet();

    public ActinCurator() {
    }

    @NotNull
    public List<ActinEntry> run(@NotNull List<ActinEntry> actinEntries) {
        List<ActinEntry> curatedActinEntry = Lists.newArrayList();

        for (ActinEntry actinEntry : actinEntries) {
            curatedActinEntry.add(curate(actinEntry));
        }
        return curatedActinEntry;
    }

    public void reportUnusedCurationEntries() {
        int unusedEntryCount = 0;
        for (CurationEntry entry : CurationFactory.MUTATION_MAPPINGS.keySet()) {
            if (!evaluatedCurationEntries.contains(entry)) {
                unusedEntryCount++;
                LOGGER.warn("Entry '{}' hasn't been used during ACTIN curation", entry);
            }
        }

        LOGGER.debug("Found {} unused ACTIN curation entries. {} keys have been requested against {} curation entries",
                unusedEntryCount,
                evaluatedCurationEntries.size(),
                CurationFactory.MUTATION_MAPPINGS.size());

        int unusedGeneCount = 0;
        for (String gene : CurationFactory.GENE_MAPPINGS.keySet()) {
            if (!evaluatedGenes.contains(gene)) {
                unusedEntryCount++;
                LOGGER.warn("Gene '{}' hasn't been used during ACTIN curation", gene);
            }
        }

        LOGGER.debug("Found {} unused ACTIN curation genes. {} genes have been requested against {} curation genes",
                unusedGeneCount,
                evaluatedGenes.size(),
                CurationFactory.GENE_MAPPINGS.size());
    }

    @NotNull
    private ActinEntry curate(@NotNull ActinEntry actinEntry) {
        List<String> parameters = actinEntry.parameters();

        String variant = Strings.EMPTY;
        String gene = Strings.EMPTY;
        if (parameters.size() == 1) {
            gene = parameters.get(0);
        } else if (parameters.size() == 2) {
            gene = parameters.get(0);
            variant = parameters.get(1);
        }

        CurationEntry entry = new CurationEntry(gene, variant);
        evaluatedCurationEntries.add(entry);
        evaluatedGenes.add(gene);

        ActinEntry curatedEntry = actinEntry;
        if (CurationFactory.MUTATION_MAPPINGS.containsKey(entry)) {
            String mappedName = CurationFactory.MUTATION_MAPPINGS.get(entry).name();
            String mappedGene = CurationFactory.MUTATION_MAPPINGS.get(entry).gene();

            LOGGER.debug("Mapping mutation '{}' on '{}' to '{}' on '{}'", entry.name(), entry.gene(), mappedName, mappedGene);
            curatedEntry = ImmutableActinEntry.builder().from(actinEntry).parameters(Lists.newArrayList(mappedGene, mappedName)).build();
        }

        if (CurationFactory.GENE_MAPPINGS.containsKey(gene)) {
            String mappedGene = CurationFactory.GENE_MAPPINGS.get(gene);
            LOGGER.debug("Mapping mutation '{}' on '{}' to '{}' on '{}'", entry.name(), entry.gene(), entry.name(), mappedGene);
            curatedEntry = ImmutableActinEntry.builder().from(actinEntry).parameters(Lists.newArrayList(mappedGene)).build();
        }

        return curatedEntry;
    }
}
