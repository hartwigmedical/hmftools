package com.hartwig.hmftools.serve.vicc.curation;

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
    private static final Set<String> NON_ONCOGENIC_INDICATORS = Sets.newHashSet("Inconclusive", "Likely Neutral");

    @NotNull
    private final Set<CurationKey> evaluatedCurationKeys = Sets.newHashSet();

    public FeatureCurator() {
    }

    @NotNull
    public List<ViccEntry> curate(@NotNull List<ViccEntry> entries) {
        List<ViccEntry> curatedViccEntries = Lists.newArrayList();

        for (ViccEntry entry : entries) {
            boolean includeEntry = potentiallyOncogenic(entry);
            List<Feature> curatedFeatures = Lists.newArrayList();
            for (Feature feature : entry.features()) {
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

    private static boolean potentiallyOncogenic(@NotNull ViccEntry entry) {
        String oncogenic = entry.association().oncogenic();
        return oncogenic == null || !NON_ONCOGENIC_INDICATORS.contains(oncogenic);
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

    @NotNull
    public Set<CurationKey> unusedCurationKeys() {
        Set<CurationKey> unusedCurationKeys = Sets.newHashSet();

        unusedCurationKeys.addAll(CurationFactory.FEATURE_BLACKLIST);
        unusedCurationKeys.addAll(CurationFactory.FEATURE_NAME_MAPPINGS.keySet());
        unusedCurationKeys.removeAll(evaluatedCurationKeys);

        return unusedCurationKeys;
    }
}
