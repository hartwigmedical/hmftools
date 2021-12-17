package com.hartwig.hmftools.serve.sources.ckb.filter;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.classification.CkbEventAndGeneExtractor;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbFilter {

    private static final Logger LOGGER = LogManager.getLogger(CkbFilter.class);

    @NotNull
    private final List<CkbFilterEntry> filters;
    @NotNull
    private final Set<CkbFilterEntry> usedFilters = Sets.newHashSet();

    public CkbFilter(@NotNull final List<CkbFilterEntry> filters) {
        this.filters = filters;
    }

    @NotNull
    public List<CkbEntry> run(@NotNull List<CkbEntry> ckbEntries) {
        List<CkbEntry> filteredCkbEntries = Lists.newArrayList();
        for (CkbEntry entry : ckbEntries) {
            if (entry.variants().size() > 1) {
                // Do not filter variants when in combination event, since this might make them a non-combined event.
                filteredCkbEntries.add(entry);
            } else if (entry.variants().isEmpty()) {
                // Always filter entries with no variants. Should never happen in practice!
                LOGGER.warn("Filtering '{}' because no variants have been defined for this entry!", entry);
            } else {
                List<Variant> filteredVariants = Lists.newArrayList();

                Variant variant = entry.variants().get(0);
                if (include(entry.type(), variant)) {
                    filteredVariants.add(variant);
                } else {
                    LOGGER.debug("Filtering variant '{}' on '{}'", variant.variant(), variant.gene().geneSymbol());
                }

                if (!filteredVariants.isEmpty()) {
                    filteredCkbEntries.add(ImmutableCkbEntry.builder().from(entry).variants(filteredVariants).build());
                }
            }
        }

        return filteredCkbEntries;
    }

    public void reportUnusedFilterEntries() {
        int unusedFilterEntryCount = 0;
        for (CkbFilterEntry entry : filters) {
            if (!usedFilters.contains(entry)) {
                unusedFilterEntryCount++;
                LOGGER.warn(" Filter entry '{}' hasn't been used for CKB filtering", entry);
            }
        }

        LOGGER.debug(" Found {} unused filter entries during CKB filtering", unusedFilterEntryCount);
    }

    private boolean include(@NotNull EventType type, @NotNull Variant variant) {
        String gene = CkbEventAndGeneExtractor.extractGene(variant);
        String event = CkbEventAndGeneExtractor.extractEvent(variant);

        for (CkbFilterEntry filterEntry : filters) {
            boolean filterMatches = isMatch(filterEntry, type, gene, event, variant.fullName());
            if (filterMatches) {
                usedFilters.add(filterEntry);
                return false;
            }
        }

        return true;
    }

    private boolean isMatch(@NotNull CkbFilterEntry filterEntry, @NotNull EventType type, @NotNull String gene,
            @NotNull String event, @NotNull String fullName) {
        switch (filterEntry.type()) {
            case FILTER_EVENT_WITH_KEYWORD: {
                return event.contains(filterEntry.value());
            }
            case FILTER_ALL_EVIDENCE_ON_GENE: {
                return gene.equals(filterEntry.value());
            }
            case FILTER_EVIDENCE_FOR_EXONS_ON_GENE: {
                return gene.equals(filterEntry.value()) && type == EventType.EXON;
            }
            case ALLOW_GENE_IN_FUSIONS_EXCLUSIVELY: {
                return gene.equals(filterEntry.value()) && type != EventType.FUSION_PAIR && type != EventType.PROMISCUOUS_FUSION;
            }
            case FILTER_SECONDARY_GENE_WHEN_FUSION_LEG: {
                return type == EventType.FUSION_PAIR && !gene.equals(filterEntry.value()) && event.contains(filterEntry.value());
            }
            case FILTER_EXACT_VARIANT_FULLNAME: {
                return fullName.equals(filterEntry.value());
            }
            default: {
                LOGGER.warn("Filter entry found with unrecognized type: {}", filterEntry);
                return false;
            }
        }
    }
}
