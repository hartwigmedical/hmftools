package com.hartwig.hmftools.serve.sources.actin.filter;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ActinFilter {

    private static final Logger LOGGER = LogManager.getLogger(ActinFilter.class);

    @NotNull
    private final List<ActinFilterEntry> filters;
    @NotNull
    private final Set<ActinFilterEntry> usedFilters = Sets.newHashSet();

    public ActinFilter(@NotNull final List<ActinFilterEntry> filters) {
        this.filters = filters;
    }

    @NotNull
    public List<ActinEntry> run(@NotNull List<ActinEntry> actinEntries) {
        List<ActinEntry> filtered = Lists.newArrayList();
        for (ActinEntry entry : actinEntries) {
            if (include(entry)) {
                filtered.add(entry);
            } else {
                LOGGER.debug("Filtering ACTIN entry '{}'", entry);
            }
        }

        return filtered;
    }

    public void reportUnusedFilterEntries() {
        int unusedFilterEntryCount = 0;
        for (ActinFilterEntry entry : filters) {
            if (!usedFilters.contains(entry)) {
                unusedFilterEntryCount++;
                LOGGER.warn(" Filter entry '{}' hasn't been used for ACTIN filtering", entry);
            }
        }

        LOGGER.debug(" Found {} unused filter entries during ACTIN filtering", unusedFilterEntryCount);
    }

    private boolean include(@NotNull ActinEntry entry) {
        for (ActinFilterEntry filterEntry : filters) {
            boolean filterMatches = isMatch(filterEntry, entry);
            if (filterMatches) {
                usedFilters.add(filterEntry);
                return false;
            }
        }

        return true;
    }

    private static boolean isMatch(@NotNull ActinFilterEntry filter, @NotNull ActinEntry entry) {
        switch (filter.type()) {
            case FILTER_EVERYTHING_FOR_RULE: {
                return ActinRule.valueOf(filter.value()) == entry.rule();
            }
            case FILTER_EVERYTHING_FOR_GENE: {
                String gene = entry.gene();
                return gene != null && gene.equals(filter.value());
            }
            case FILTER_VARIANT_ON_GENE: {
                String gene = entry.gene();
                if (gene == null) {
                    return false;
                } else {
                    String evaluation = gene + " " + entry.mutation();
                    return evaluation.equals(filter.value());
                }
            } default: {
                LOGGER.warn("Filter entry found with unrecognized type: {}", filter);
                return false;
            }
        }
    }
}
