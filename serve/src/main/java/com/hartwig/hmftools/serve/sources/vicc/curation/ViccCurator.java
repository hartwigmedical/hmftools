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

public class ViccCurator {

    private static final Logger LOGGER = LogManager.getLogger(ViccCurator.class);

    @NotNull
    private final Set<CurationKey> evaluatedCurationKeys = Sets.newHashSet();

    public ViccCurator() {
    }

    @NotNull
    public List<ViccEntry> run(@NotNull List<ViccEntry> entries) {
        List<ViccEntry> curatedViccEntries = Lists.newArrayList();

        for (ViccEntry entry : entries) {
            List<Feature> curatedFeatures = Lists.newArrayList();
            boolean includeEntry = true;
            for (Feature feature : entry.features()) {
                if (feature.name().startsWith("PAX8-PPAR")) {
                    int x = 1;
                }
                Feature curatedFeature = curate(entry, feature);
                if (curatedFeature != null) {
                    curatedFeatures.add(curatedFeature);
                } else {
                    includeEntry = false;
                }
            }

            if (includeEntry) {
                curatedViccEntries.add(ImmutableViccEntry.builder().from(entry).features(curatedFeatures).build());
            }
        }
        return curatedViccEntries;
    }

    public void reportUnusedCurationEntries() {
        int unusedKeyCount = 0;
        for (CurationKey key : CurationFactory.FEATURE_BLACKLIST) {
            if (!evaluatedCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.warn("Key '{}' hasn't been used during VICC blacklist curation", key);
            }
        }

        for (CurationKey key : CurationFactory.FEATURE_NAME_MAPPINGS.keySet()) {
            if (!evaluatedCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.warn("Key '{}' hasn't been used during VICC name mapping curation", key);
            }
        }

        int totalKeyCount = CurationFactory.FEATURE_BLACKLIST.size() + CurationFactory.FEATURE_NAME_MAPPINGS.keySet().size();
        LOGGER.debug("Found {} unused VICC curation entries. {} keys have been requested against {} blacklist entries",
                unusedKeyCount,
                evaluatedCurationKeys.size(),
                totalKeyCount);
    }

    @VisibleForTesting
    @Nullable
    Feature curate(@NotNull ViccEntry entry, @NotNull Feature feature) {
        CurationKey key = new CurationKey(entry.source(), feature.geneSymbol(), entry.transcriptId(), feature.name());
        evaluatedCurationKeys.add(key);

        if (CurationFactory.FEATURE_BLACKLIST.contains(key)) {
            LOGGER.debug("Blacklisting feature '{}' for gene {} in {}", feature.name(), feature.geneSymbol(), entry.source());
            return null;
        } else {
            String mappedFeatureName = CurationFactory.FEATURE_NAME_MAPPINGS.get(key);
            if (mappedFeatureName != null) {
                LOGGER.debug("Mapping feature '{}' to '{}' for gene {} in {}",
                        feature.name(),
                        mappedFeatureName,
                        feature.geneSymbol(),
                        entry.source());
                return ImmutableFeature.builder().from(feature).name(mappedFeatureName).build();
            }
        }

        return feature;
    }
}
