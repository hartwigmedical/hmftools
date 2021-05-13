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
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbFilter {

    private static final Logger LOGGER = LogManager.getLogger(CkbFilter.class);

    @NotNull
    private final Set<String> usedFilterKeywords = Sets.newHashSet();
    @NotNull
    private final Set<String> usedAllFilterGenes = Sets.newHashSet();
    @NotNull
    private final Set<String> usedExonFilterGenes = Sets.newHashSet();
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
        int unusedKeywordCount = 0;
        for (String keyword : FilterFactory.VARIANT_KEYWORDS_TO_FILTER) {
            if (!usedFilterKeywords.contains(keyword)) {
                unusedKeywordCount++;
                LOGGER.warn(" Keyword '{}' hasn't been used for CKB filtering", keyword);
            }
        }

        LOGGER.debug(" Found {} unused keywords during CKB filtering", unusedKeywordCount);

        int unusedAllGeneCount = 0;
        for (String gene : FilterFactory.GENES_FOR_WHICH_TO_FILTER_ALL) {
            if (!usedAllFilterGenes.contains(gene)) {
                unusedAllGeneCount++;
                LOGGER.warn(" Gene '{}' hasn't been used for CKB all gene event filtering", gene);
            }
        }

        LOGGER.debug(" Found {} unused genes during CKB all gene event filtering", unusedAllGeneCount);

        int unusedExonGeneCount = 0;
        for (String gene : FilterFactory.GENES_FOR_WHICH_TO_FILTER_EXON_EVENTS) {
            if (!usedExonFilterGenes.contains(gene)) {
                unusedExonGeneCount++;
                LOGGER.warn(" Gene '{}' hasn't been used for exon gene event CKB filtering", gene);
            }
        }

        LOGGER.debug(" Found {} unused genes during exon gene event CKB filtering", unusedExonGeneCount);
    }

    private boolean include(@NotNull EventType type, @NotNull Variant variant) {
        String event = ckbEventAndGeneExtractor.extractEvent(variant);
        for (String keywordToFilter : FilterFactory.VARIANT_KEYWORDS_TO_FILTER) {
            if (event.contains(keywordToFilter)) {
                usedFilterKeywords.add(keywordToFilter);
                return false;
            }
        }

        String gene = ckbEventAndGeneExtractor.extractGene(variant);
        if (FilterFactory.GENES_FOR_WHICH_TO_FILTER_ALL.contains(gene)) {
            usedAllFilterGenes.add(gene);
            return false;
        }

        if (type == EventType.EXON) {
            if (FilterFactory.GENES_FOR_WHICH_TO_FILTER_EXON_EVENTS.contains(gene)) {
                usedExonFilterGenes.add(gene);
                return false;
            }
        }

        // We don't want to include evidence on genes that are unmappable between ref genome versions.
        if (CkbConstants.UNMAPPABLE_GENES.contains(gene)) {
            return false;
        }

        // For some genes we only accept evidence in case the evidence is on a fusion.
        if (type != EventType.FUSION_PAIR && type != EventType.PROMISCUOUS_FUSION) {
            if (CkbConstants.EXCLUSIVE_FUSION_GENES.contains(gene)) {
                return false;
            }
        }

        // In case a leg of a fusion does not exist in our exome definition, or mapping, we have to remove it.
        if (type == EventType.FUSION_PAIR) {
            for (String fusionLeg : CkbConstants.UNRESOLVABLE_FUSION_LEGS) {
                if (event.contains(fusionLeg) && !gene.equals(fusionLeg)) {
                    return false;
                }
            }

            for (String unmappableGene : CkbConstants.UNMAPPABLE_GENES) {
                if (event.contains(unmappableGene) && !gene.equals(unmappableGene)) {
                    return false;
                }
            }
        }

        return true;
    }
}
