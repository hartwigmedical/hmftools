package com.hartwig.hmftools.serve.vicc.range;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneLevelEventExtractor {
    private static final Logger LOGGER = LogManager.getLogger(GeneLevelEventExtractor.class);

    private static final Set<String> GENE = Sets.newHashSet("mut", "mutant");

    private static final Set<String> GENE_ACTIVATION = Sets.newHashSet("Gain-of-function Mutations", "act mut");

    private static final Set<String> GENE_INACTIVATION = Sets.newHashSet("Truncating Mutations", "inact mut", "loss");

    @NotNull
    public Map<Feature, String> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> geneLevelEventsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {

            if (GENE_ACTIVATION.contains(feature.name()) || GENE_ACTIVATION.contains(feature.biomarkerType())) {
                geneLevelEventsPerFeature.put(feature, "gain of " + feature.geneSymbol());
            } else if (GENE_INACTIVATION.contains(feature.name()) || GENE_INACTIVATION.contains(feature.biomarkerType())) {
                geneLevelEventsPerFeature.put(feature, "loss of " + feature.geneSymbol());
            } else if (GENE.contains(feature.biomarkerType())) { //TODO: determine gain of loss function
                geneLevelEventsPerFeature.put(feature, " " + feature.geneSymbol());
            }
        }

        return geneLevelEventsPerFeature;
    }
}
