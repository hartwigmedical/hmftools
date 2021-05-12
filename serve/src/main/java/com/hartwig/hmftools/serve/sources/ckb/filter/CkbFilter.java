package com.hartwig.hmftools.serve.sources.ckb.filter;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.classification.CkbConstants;
import com.hartwig.hmftools.ckb.classification.CkbEventAndGeneExtractor;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbFilter {

    private static final Logger LOGGER = LogManager.getLogger(CkbFilter.class);

    @NotNull
    private final Set<String> usedFilterKeywords = Sets.newHashSet();
    @NotNull
    private final Set<String> usedFilterGenes = Sets.newHashSet();
    @NotNull
    private final CkbEventAndGeneExtractor ckbEventAndGeneExtractor = new CkbEventAndGeneExtractor();

    public CkbFilter() {
    }

    @NotNull
    public List<CkbEntry> run(@NotNull List<CkbEntry> ckbEntries) {
        List<CkbEntry> filteredCkbEntries = Lists.newArrayList();
        for (CkbEntry entry : ckbEntries) {
            if (entry.variants().size() > 1) {
                // Do not filter variants when in combination event, since this might make them a non-combined event.
                filteredCkbEntries.add(entry);
            } else if (entry.variants().size() == 0) {
                // Always filter entries with no variants. Should never happen in practice!
                LOGGER.warn("Filtering '{}' because no variants have been defined for this entry!", entry);
            } else {
                List<Variant> filteredVariants = Lists.newArrayList();

                Variant variant = entry.variants().get(0);
                if (include(variant)) {
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
        int unusedKeywordCount = 0;
        for (String keyword : FilterFactory.VARIANT_KEYWORDS_TO_FILTER) {
            if (!usedFilterKeywords.contains(keyword)) {
                unusedKeywordCount++;
                LOGGER.warn(" Keyword '{}' hasn't been used for CKB filtering", keyword);
            }
        }

        LOGGER.debug(" Found {} unused keywords during CKB filtering", unusedKeywordCount);

        int unusedGeneCount = 0;
        for (String gene : FilterFactory.GENES_TO_FILTER) {
            if (!usedFilterGenes.contains(gene)) {
                unusedGeneCount++;
                LOGGER.warn(" Gene '{}' hasn't been used for CKB filtering", gene);
            }
        }

        LOGGER.debug(" Found {} unused genes during CKB filtering", unusedGeneCount);
    }

    private boolean include(@NotNull Variant variant) {
        String event = ckbEventAndGeneExtractor.extractEvent(variant);
        for (String keywordToFilter : FilterFactory.VARIANT_KEYWORDS_TO_FILTER) {
            if (event.contains(keywordToFilter)) {
                usedFilterKeywords.add(keywordToFilter);
                return false;
            }
        }

        String gene = ckbEventAndGeneExtractor.extractGene(variant);
        if (FilterFactory.GENES_TO_FILTER.contains(gene)) {
            usedFilterGenes.add(gene);
            return false;
        }

        // We don't want to include evidence on genes that are unmappable between ref genome versions.
        if (CkbConstants.UNMAPPABLE_GENES.contains(gene)) {
            return false;
        }

        return true;
    }
}
