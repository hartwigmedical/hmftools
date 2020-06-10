package com.hartwig.hmftools.serve.vicc.curation;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class FeatureCurator {

    private static final Logger LOGGER = LogManager.getLogger(FeatureCurator.class);

    @VisibleForTesting
    static final Map<String, List<FeatureNameMapping>> ONCOKB_FEATURE_NAME_MAPPINGS_PER_GENE = Maps.newHashMap();

    static {
        ONCOKB_FEATURE_NAME_MAPPINGS_PER_GENE.put("EPAS1",
                Lists.newArrayList(ImmutableFeatureNameMapping.builder()
                        .originalFeatureName("533_534del")
                        .curatedFeatureName("I533_P534del")
                        .build()));
    }

    private FeatureCurator() {
    }

    @NotNull
    public static Feature curate(@NotNull ViccSource source, @NotNull Feature feature) {
        if (source == ViccSource.ONCOKB) {
            List<FeatureNameMapping> mappings = ONCOKB_FEATURE_NAME_MAPPINGS_PER_GENE.get(feature.geneSymbol());
            if (mappings != null) {
                for (FeatureNameMapping mapping : mappings) {
                    if (mapping.originalFeatureName().equals(feature.name())) {
                        LOGGER.debug("Curating feature '{}' to '{}' for gene {} in {}",
                                mapping.originalFeatureName(),
                                mapping.curatedFeatureName(),
                                feature.geneSymbol(),
                                source);
                        return ImmutableFeature.builder().from(feature).name(mapping.curatedFeatureName()).build();
                    }
                }
            }
        }

        return feature;
    }
}
