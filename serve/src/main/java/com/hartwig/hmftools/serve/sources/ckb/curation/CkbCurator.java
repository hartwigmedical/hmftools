package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableGene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbCurator {

    private static final Logger LOGGER = LogManager.getLogger(CkbCurator.class);

    @NotNull
    private final Set<String> evaluatedGenes = Sets.newHashSet();
    @NotNull
    private final Set<CurationEntry> evaluatedCurationEntries = Sets.newHashSet();

    public CkbCurator() {

    }

    @NotNull
    public List<CkbEntry> run(@NotNull List<CkbEntry> ckbEntries) {
        List<CkbEntry> curatedCkbEntries = Lists.newArrayList();

        for (CkbEntry ckbEntry : ckbEntries) {
            List<Variant> curatedVariants = Lists.newArrayList();
            if (ckbEntry.variants().size() == 1) {
                curatedVariants.add(curate(ckbEntry.variants()));

                curatedCkbEntries.add(ImmutableCkbEntry.builder().from(ckbEntry).variants(curatedVariants).build());

            }
        }

        return curatedCkbEntries;
    }

    public void reportUnusedCurationEntries() {
        int unusedEntryCount = 0;
        for (CurationEntry entry : CurationFactory.MUTATION_MAPPINGS.keySet()) {
            if (!evaluatedCurationEntries.contains(entry)) {
                unusedEntryCount++;
                LOGGER.warn("Entry '{}' hasn't been used during CKB curation", entry);
            }
        }

        LOGGER.info("Found {} unused CKB curation entries. {} keys have been requested against {} curation entries",
                unusedEntryCount,
                evaluatedCurationEntries.size(),
                CurationFactory.MUTATION_MAPPINGS.size());

        int unusedGeneCount = 0;
        for (String gene : CurationFactory.GENE_MAPPINGS.keySet()) {
            if (!evaluatedGenes.contains(gene)) {
                unusedEntryCount++;
                LOGGER.warn("Gene '{}' hasn't been used during CKB curation", gene);
            }
        }

        LOGGER.info("Found {} unused CKB curation genes. {} genes have been requested against {} curation genes",
                unusedGeneCount,
                evaluatedGenes.size(),
                CurationFactory.GENE_MAPPINGS.size());
    }

    @NotNull
    private Variant curate(@NotNull List<Variant> variants) {
        Variant variant = variants.get(0);
        CurationEntry entry = new CurationEntry(variant.gene().geneSymbol(), variant.variant());
        evaluatedCurationEntries.add(entry);

        String gene = variant.gene().geneSymbol();
        evaluatedGenes.add(gene);

        Variant curatedMutation = variant;
        if (CurationFactory.MUTATION_MAPPINGS.containsKey(entry)) {
            String mappedName = CurationFactory.MUTATION_MAPPINGS.get(entry).name();
            String mappedGene = CurationFactory.MUTATION_MAPPINGS.get(entry).gene();

            LOGGER.debug("Mapping mutation '{}' on '{}' to '{}' on '{}'", entry.name(), entry.gene(), mappedName, mappedGene);
            curatedMutation = ImmutableVariant.builder()
                    .from(curatedMutation)
                    .gene(ImmutableGene.builder().from(variant.gene()).geneSymbol(mappedGene).build())
                    .variant(mappedName)
                    .build();
        }

        if (CurationFactory.GENE_MAPPINGS.containsKey(gene)) {
            String mappedGene = CurationFactory.GENE_MAPPINGS.get(gene);
            LOGGER.debug("Mapping mutation '{}' on '{}' to '{}' on '{}'", entry.name(), entry.gene(), entry.name(), mappedGene);
            curatedMutation = ImmutableVariant.builder()
                    .from(curatedMutation)
                    .gene(ImmutableGene.builder().from(variant.gene()).geneSymbol(mappedGene).build())
                    .build();
        }

        return curatedMutation;
    }
}
