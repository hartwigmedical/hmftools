package com.hartwig.hmftools.vicc.annotation;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class MutationTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(MutationTypeExtractor.class);

    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(new ProteinAnnotationExtractor());

    @NotNull
    public static MutationType extractType(@NotNull Feature feature) {
        String gene = feature.geneSymbol();
        if (gene == null) {
            LOGGER.debug("Skipping extraction for '{}' since gene is missing", feature.name());
            return MutationType.UNKNOWN;
        } else {
            return CLASSIFIER.determineType(gene, feature.name());
        }
    }
}
