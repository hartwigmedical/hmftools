package com.hartwig.hmftools.serve.vicc.curation;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FeatureCurator {

    private static final Logger LOGGER = LogManager.getLogger(FeatureCurator.class);

    @VisibleForTesting
    static final Map<CurationKey, List<FeatureNameMapping>> ONCOKB_FEATURE_NAME_MAPPINGS = Maps.newHashMap();

    static {
        ONCOKB_FEATURE_NAME_MAPPINGS.put(new CurationKey("EPAS1", "ENST00000263734"),
                Lists.newArrayList(ImmutableFeatureNameMapping.builder()
                        .originalFeatureName("533_534del")
                        .curatedFeatureName("I533_P534del")
                        .build()));

        ONCOKB_FEATURE_NAME_MAPPINGS.put(new CurationKey("KIT", "ENST00000288135"),
                Lists.newArrayList(ImmutableFeatureNameMapping.builder()
                        .originalFeatureName("V559del")
                        .curatedFeatureName("V560del")
                        .build()));

        ONCOKB_FEATURE_NAME_MAPPINGS.put(new CurationKey("PTEN", "ENST00000371953"),
                Lists.newArrayList(ImmutableFeatureNameMapping.builder()
                        .originalFeatureName("I32del")
                        .curatedFeatureName("I33del")
                        .build()));
    }

    private FeatureCurator() {
    }

    @NotNull
    public static Feature curate(@NotNull ViccEntry entry, @NotNull Feature feature) {
        CurationKey key = new CurationKey(feature.geneSymbol(), entry.transcriptId());
        if (entry.source() == ViccSource.ONCOKB) {
            List<FeatureNameMapping> mappings = ONCOKB_FEATURE_NAME_MAPPINGS.get(key);
            if (mappings != null) {
                for (FeatureNameMapping mapping : mappings) {
                    if (mapping.originalFeatureName().equals(feature.name())) {
                        LOGGER.debug("Curating feature '{}' to '{}' for gene {} in {}",
                                mapping.originalFeatureName(),
                                mapping.curatedFeatureName(),
                                feature.geneSymbol(),
                                entry);
                        return ImmutableFeature.builder().from(feature).name(mapping.curatedFeatureName()).build();
                    }
                }
            }
        }

        return feature;
    }

    @VisibleForTesting
    static class CurationKey {
        @NotNull
        private final String gene;
        @Nullable
        private final String transcript;

        public CurationKey(@NotNull final String gene, @Nullable final String transcript) {
            this.gene = gene;
            this.transcript = transcript;
        }

        @VisibleForTesting
        @NotNull
        String gene() {
            return gene;
        }

        @VisibleForTesting
        @Nullable
        String transcript() {
            return transcript;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final CurationKey that = (CurationKey) o;
            return gene.equals(that.gene) && Objects.equals(transcript, that.transcript);
        }

        @Override
        public int hashCode() {
            return Objects.hash(gene, transcript);
        }
    }
}
