package com.hartwig.hmftools.serve.sources.vicc.filter;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViccFilter {

    private static final Logger LOGGER = LogManager.getLogger(ViccFilter.class);

    private static final Set<String> NON_ONCOGENIC_INDICATORS = Sets.newHashSet("Inconclusive", "Likely Neutral");

    @NotNull
    private final Set<String> filteredKeywords = Sets.newHashSet();
    @NotNull
    private final Set<String> filteredFeatures = Sets.newHashSet();
    @NotNull
    private final Set<FilterKey> filteredKeys = Sets.newHashSet();

    public ViccFilter() {
    }

    @NotNull
    public List<ViccEntry> run(@NotNull List<ViccEntry> entries) {
        List<ViccEntry> filteredViccEntries = Lists.newArrayList();

        for (ViccEntry entry : entries) {
            // We exclude any non-oncogenic event for further processing.
            if (potentiallyOncogenic(entry)) {
                List<Feature> filteredFeatures = Lists.newArrayList();
                boolean hasFilteredFeatures = false;
                for (Feature feature : entry.features()) {
                    if (include(entry.source(), feature)) {
                        filteredFeatures.add(feature);
                    } else {
                        hasFilteredFeatures = true;
                      //  LOGGER.debug("Filtering feature '{}' on '{}'", feature.name(), feature.geneSymbol());
                    }
                }

                // Only exclude VICC entries with no features in case at least one feature has been filtered.
                if (!filteredFeatures.isEmpty() || !hasFilteredFeatures) {
                    filteredViccEntries.add(ImmutableViccEntry.builder().from(entry).features(filteredFeatures).build());
                }
            }
        }

        return filteredViccEntries;
    }

    public void reportUnusedFilterEntries() {
        int unusedKeywordCount = 0;
        for (String keyword : FilterFactory.FEATURE_KEYWORDS_TO_FILTER) {
            if (!filteredKeywords.contains(keyword)) {
                unusedKeywordCount++;
                LOGGER.warn("Keyword '{}' hasn't been used for VICC filtering", keyword);
            }
        }

        int unusedFeatureCount = 0;
        for (String feature : FilterFactory.FEATURES_TO_FILTER) {
            if (!filteredFeatures.contains(feature)) {
                unusedFeatureCount++;
                LOGGER.warn("Feature '{}' hasn't been used for VICC filtering", feature);
            }
        }

        int unusedKeyCount = 0;
        for (FilterKey key : FilterFactory.FEATURE_KEYS_TO_FILTER) {
            if (!filteredKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.warn("Key '{}' hasn't been used for VICC filtering", key);
            }
        }

        LOGGER.debug("Found {} unused keywords, {} unused features and {} unused keys during VICC filtering",
                unusedKeywordCount,
                unusedFeatureCount,
                unusedKeyCount);
    }

    @VisibleForTesting
    boolean include(@NotNull ViccSource source, @NotNull Feature feature) {
        String featureName = feature.name();
        for (String keywordToFilter : FilterFactory.FEATURE_KEYWORDS_TO_FILTER) {
            if (featureName.contains(keywordToFilter)) {
                filteredKeywords.add(keywordToFilter);
                return false;
            }
        }

        for (String featureToFilter : FilterFactory.FEATURES_TO_FILTER) {
            if (featureName.equals(featureToFilter)) {
                filteredFeatures.add(featureToFilter);
                return false;
            }
        }

        String gene = feature.geneSymbol();
        if (gene != null) {
            FilterKey filterKey = new FilterKey(source, gene, featureName);
            if (FilterFactory.FEATURE_KEYS_TO_FILTER.contains(filterKey)) {
                filteredKeys.add(filterKey);
                return false;
            }
        }

        return true;
    }

    private static boolean potentiallyOncogenic(@NotNull ViccEntry entry) {
        String oncogenic = entry.association().oncogenic();
        return oncogenic == null || !NON_ONCOGENIC_INDICATORS.contains(oncogenic);
    }
}
