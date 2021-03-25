package com.hartwig.hmftools.serve.sources.ckb.filter;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CKBFilter {

    private static final Logger LOGGER = LogManager.getLogger(CKBFilter.class);

    @NotNull
    private final Set<String> filteredMutations = Sets.newHashSet();


    public CKBFilter() {
    }


    @NotNull
    public List<CkbEntry> run(@NotNull List<CkbEntry> entries) {

        List<CkbEntry> filteredTrials = Lists.newArrayList();
        for (CkbEntry entry : entries) {
            if (entry.variants().size() == 1) {
                List<Variant> filteredVariants = Lists.newArrayList();

                if (include(entry.variants().get(0))) {
                    filteredVariants.add(entry.variants().get(0));
                } else {
                    LOGGER.info("Filtering feature '{}' on '{}'", entry.variants().get(0).variant(), entry.variants().get(0).gene().geneSymbol());
                }

                if (!filteredVariants.isEmpty()) {
                    filteredTrials.add(ImmutableCkbEntry.builder().from(entry).variants(filteredVariants).build());
                }

            }

        }


        return filteredTrials;

    }

    public void reportUnusedFilterEntries() {
        int unusedKeywordCount = 0;
        for (String keyword : FilterFactory.PROFILE_NAME_TO_FILTER) {
            if (!filteredMutations.contains(keyword)) {
                unusedKeywordCount++;
                LOGGER.debug("Keyword '{}' hasn't been used for CKB filtering", keyword);
            }
        }

               LOGGER.debug("Found {} unused keywords during CKB filtering",
                unusedKeywordCount);
    }

    @VisibleForTesting
    boolean include(@NotNull Variant variant) {
        String profileName = variant.variant();
        for (String keywordToFilter : FilterFactory.PROFILE_NAME_TO_FILTER) {
            if (profileName.contains(keywordToFilter)) {
                filteredMutations.add(keywordToFilter);
                return false;
            }
        }

        return true;
    }
}
