package com.hartwig.hmftools.serve.sources.actin.filter;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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
        List<ActinEntry> filteredActinEntries = Lists.newArrayList();
        for (ActinEntry entry : actinEntries) {
            List<String> filteredVariants = Lists.newArrayList();

            String variant = Strings.EMPTY;
            String gene = Strings.EMPTY;
            if (entry.parameters().size() == 2) {
                gene = entry.parameters().get(0);
                variant = entry.parameters().get(1);
            }

            if (entry.parameters().size() == 1) {
                gene = entry.parameters().get(0);
            }

            for (EventType eventType : entry.type()) {
                if (include(eventType, entry)) {
                    filteredVariants.add(variant);
                } else {
                    LOGGER.debug("Filtering variant '{}' on '{}'", variant, gene);
                }

                if (!filteredVariants.isEmpty()) {
                    filteredActinEntries.add(ImmutableActinEntry.builder().from(entry).parameters(entry.parameters()).build());
                }
            }

        }

        return filteredActinEntries;
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

    private boolean include(@NotNull EventType type, @NotNull ActinEntry entry) {
        for (ActinFilterEntry filterEntry : filters) {
            boolean filterMatches = isMatch(filterEntry, entry);
            if (filterMatches) {
                usedFilters.add(filterEntry);
                return false;
            }
        }

        return true;
    }

    private boolean isMatch(@NotNull ActinFilterEntry filterEntry, @NotNull ActinEntry entry) {

        String variant = Strings.EMPTY;
        String gene = Strings.EMPTY;

        if (entry.parameters().size() == 2) {
            gene = entry.parameters().get(0);
            variant = entry.parameters().get(1);
        }

        if (entry.parameters().size() == 1) {
            gene = entry.parameters().get(0);
        }

        String combined = gene + ", " + variant;

        switch (filterEntry.type()) {
            case FILTER_EXACT_VARIANT_FULLNAME: {
                return combined.equals(filterEntry.value());
            }
            default: {
                LOGGER.warn("Filter entry found with unrecognized type: {}", filterEntry);
                return false;
            }
        }
    }

}
