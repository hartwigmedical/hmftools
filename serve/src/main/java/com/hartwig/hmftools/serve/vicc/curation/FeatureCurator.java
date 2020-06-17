package com.hartwig.hmftools.serve.vicc.curation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FeatureCurator {

    private static final Logger LOGGER = LogManager.getLogger(FeatureCurator.class);

    @NotNull
    private final Map<ViccSource, Set<CurationKey>> evaluatedCurationKeysPerSource = Maps.newHashMap();

    public FeatureCurator() {
    }

    @Nullable
    public Feature curate(@NotNull ViccEntry entry, @NotNull Feature feature) {
        CurationKey key = new CurationKey(feature.geneSymbol(), entry.transcriptId(), feature.name());
        Set<CurationKey> keys = evaluatedCurationKeysPerSource.get(entry.source());
        if (keys == null) {
            keys = Sets.newHashSet();
        }
        keys.add(key);
        evaluatedCurationKeysPerSource.put(entry.source(), keys);

        Set<CurationKey> blacklistForSource = blacklistForSource(entry.source());
        if (blacklistForSource.contains(key)) {
            LOGGER.debug("Blacklisting feature '{}' for gene {} in {}", feature.name(), feature.geneSymbol(), entry.source());
            return null;
        } else {
            Map<CurationKey, String> mappingsForSource = mappingsForSource(entry.source());
            String mappedFeatureName = mappingsForSource.get(key);
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
    private static Set<CurationKey> blacklistForSource(@NotNull ViccSource source) {
        switch (source) {
            case CIVIC:
                return CurationFactory.CIVIC_FEATURE_BLACKLIST;
            case JAX:
                return CurationFactory.JAX_FEATURE_BLACKLIST;
            case ONCOKB:
                return CurationFactory.ONCOKB_FEATURE_BLACKLIST;
            default:
                return Sets.newHashSet();
        }
    }

    @NotNull
    private static Map<CurationKey, String> mappingsForSource(@NotNull ViccSource source) {
        switch (source) {
            case CIVIC:
                return CurationFactory.CIVIC_FEATURE_NAME_MAPPINGS;
            case JAX:
                return CurationFactory.JAX_FEATURE_NAME_MAPPINGS;
            case ONCOKB:
                return CurationFactory.ONCOKB_FEATURE_NAME_MAPPINGS;
            default:
                return Maps.newHashMap();
        }
    }

    @NotNull
    public Map<ViccSource, Set<CurationKey>> unusedCurationKeysPerSource() {
        Map<ViccSource, Set<CurationKey>> unusedCurationKeysPerSource = Maps.newHashMap();
        for (Map.Entry<ViccSource, Set<CurationKey>> entry : evaluatedCurationKeysPerSource.entrySet()) {
            ViccSource source = entry.getKey();
            Set<CurationKey> unusedKeys = Sets.newHashSet();
            unusedKeys.addAll(blacklistForSource(source));
            unusedKeys.addAll(mappingsForSource(source).keySet());
            unusedKeys.removeAll(entry.getValue());
            unusedCurationKeysPerSource.put(source, unusedKeys);
        }

        return unusedCurationKeysPerSource;
    }
}
