package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FeatureCurator {

    private static final Logger LOGGER = LogManager.getLogger(FeatureCurator.class);

    @NotNull
    private final Set<FeatureCurationKey> evaluatedFeatureCurationKeys = Sets.newHashSet();

    public FeatureCurator() {
    }

    @NotNull
    public List<ViccEntry> run(@NotNull List<ViccEntry> entries) {
        List<ViccEntry> curatedViccEntries = Lists.newArrayList();

        for (ViccEntry entry : entries) {
            List<Feature> curatedFeatures = Lists.newArrayList();
            boolean includeEntry = true;
            for (Feature feature : entry.features()) {
                Feature curatedFeature = curate(entry, feature);
                if (curatedFeature != null) {
                    curatedFeatures.add(curatedFeature);
                } else {
                    includeEntry = false;
                }
            }

            // We want to exclude entries completely with at least one blacklisted feature.
            if (includeEntry) {
                curatedViccEntries.add(ImmutableViccEntry.builder().from(entry).features(curatedFeatures).build());
            }
        }
        return curatedViccEntries;
    }

    public void reportUnusedCurationKeys() {
        int unusedKeyCount = 0;
        for (FeatureCurationKey key : FeatureCurationFactory.FEATURE_BLACKLIST) {
            if (!evaluatedFeatureCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.debug("Key '{}' hasn't been used during VICC blacklist curation", key);
            }
        }

        for (FeatureCurationKey key : FeatureCurationFactory.FEATURE_MAPPINGS.keySet()) {
            if (!evaluatedFeatureCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.debug("Key '{}' hasn't been used during VICC name mapping curation", key);
            }
        }

        int totalKeyCount = FeatureCurationFactory.FEATURE_BLACKLIST.size() + FeatureCurationFactory.FEATURE_MAPPINGS.keySet().size();
        LOGGER.debug("Found {} unused VICC feature curation entries. {} keys have been requested against {} entries",
                unusedKeyCount,
                evaluatedFeatureCurationKeys.size(),
                totalKeyCount);
    }

    @VisibleForTesting
    @Nullable
    Feature curate(@NotNull ViccEntry entry, @NotNull Feature feature) {
        FeatureCurationKey key = new FeatureCurationKey(entry.source(), feature.geneSymbol(), entry.transcriptId(), feature.name());
        evaluatedFeatureCurationKeys.add(key);

        if (FeatureCurationFactory.FEATURE_BLACKLIST.contains(key)) {
            LOGGER.debug("Blacklisting feature '{}' for gene {} in {}", feature.name(), feature.geneSymbol(), entry.source());
            return null;
        } else if (FeatureCurationFactory.FEATURE_MAPPINGS.containsKey(key)) {
            String mappedGeneSymbol = FeatureCurationFactory.FEATURE_MAPPINGS.get(key).geneSymbol();
            String mappedFeatureName = FeatureCurationFactory.FEATURE_MAPPINGS.get(key).featureName();

            LOGGER.debug("Mapping feature '{}' to '{}' for gene {} in {}",
                    feature.name(),
                    mappedFeatureName,
                    feature.geneSymbol(),
                    entry.source());
            return ImmutableFeature.builder().from(feature).geneSymbol(mappedGeneSymbol).name(mappedFeatureName).build();
        }

        return feature;
    }
}
