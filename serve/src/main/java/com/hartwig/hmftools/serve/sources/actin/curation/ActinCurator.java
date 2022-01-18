package com.hartwig.hmftools.serve.sources.actin.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ActinCurator {

    private static final Logger LOGGER = LogManager.getLogger(ActinCurator.class);

    @NotNull
    private final Set<String> evaluatedGenes = Sets.newHashSet();

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
        int unusedGeneCount = 0;
        for (String gene : CurationFactory.GENE_MAPPINGS.keySet()) {
            if (!evaluatedGenes.contains(gene)) {
                unusedGeneCount++;
                LOGGER.warn("Gene '{}' hasn't been used during ACTIN curation", gene);
            }
        }

        LOGGER.debug("Found {} unused ACTIN curation genes. {} genes have been requested against {} curation genes",
                unusedGeneCount,
                evaluatedGenes.size(),
                CurationFactory.GENE_MAPPINGS.size());
    }

    @NotNull
    private ActinEntry curate(@NotNull ActinEntry entry) {
        evaluatedGenes.add(entry.gene());

        if (!CurationFactory.GENE_MAPPINGS.containsKey(entry.gene())) {
            return entry;
        }

        String mappedGene = CurationFactory.GENE_MAPPINGS.get(entry.gene());
        LOGGER.debug("Mapping gene '{}' to '{}'", entry.gene(), mappedGene);

        return ImmutableActinEntry.builder().from(entry).gene(mappedGene).build();
    }
}
