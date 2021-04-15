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
    private final Set<CurationEntry> evaluatedCurationEntries = Sets.newHashSet();

    public CkbCurator() {
    }

    @NotNull
    public List<CkbEntry> run(@NotNull List<CkbEntry> ckbEntries) {
        List<CkbEntry> curatedCkbEntries = Lists.newArrayList();

        for (CkbEntry ckbEntry : ckbEntries) {
            List<Variant> curatedVariants = Lists.newArrayList();
            for (Variant variant : ckbEntry.variants()) {
                curatedVariants.add(curate(variant));
            }
            curatedCkbEntries.add(ImmutableCkbEntry.builder().from(ckbEntry).variants(curatedVariants).build());
        }

        return curatedCkbEntries;
    }

    public void reportUnusedCurationEntries() {
        int unusedEntryCount = 0;
        for (CurationEntry entry : CurationFactory.VARIANT_MAPPINGS.keySet()) {
            if (!evaluatedCurationEntries.contains(entry)) {
                unusedEntryCount++;
                LOGGER.warn(" Entry '{}' hasn't been used during CKB curation", entry);
            }
        }

        LOGGER.debug(" Found {} unused CKB curation entries. {} keys have been requested against {} curation entries",
                unusedEntryCount,
                evaluatedCurationEntries.size(),
                CurationFactory.VARIANT_MAPPINGS.size());
    }

    @NotNull
    private Variant curate(@NotNull Variant variant) {
        CurationEntry entry = new CurationEntry(variant.gene().geneSymbol(), variant.variant());
        evaluatedCurationEntries.add(entry);

        Variant curatedMutation = variant;
        if (CurationFactory.VARIANT_MAPPINGS.containsKey(entry)) {
            String mappedVariant = CurationFactory.VARIANT_MAPPINGS.get(entry).variant();
            String mappedGeneSymbol = CurationFactory.VARIANT_MAPPINGS.get(entry).geneSymbol();

            LOGGER.debug("Mapping mutation '{}' on '{}' to '{}' on '{}'",
                    entry.variant(),
                    entry.geneSymbol(),
                    mappedVariant,
                    mappedGeneSymbol);

            curatedMutation = ImmutableVariant.builder()
                    .from(curatedMutation)
                    .gene(ImmutableGene.builder().from(variant.gene()).geneSymbol(mappedGeneSymbol).build())
                    .variant(mappedVariant)
                    .build();
        }

        return curatedMutation;
    }
}
