package com.hartwig.hmftools.vicc.annotation;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.MutationType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class MutationTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(MutationTypeExtractor.class);

    @NotNull
    private final EventClassifier classifier;

    @NotNull
    public static MutationTypeExtractor buildProductionExtractor() {
        EventClassifier eventClassifier = EventClassifierFactory.buildClassifier(new ProteinAnnotationExtractor());
        return new MutationTypeExtractor(eventClassifier);
    }

    private MutationTypeExtractor(@NotNull final EventClassifier classifier) {
        this.classifier = classifier;
    }

    @NotNull
    public MutationType extractType(@Nullable String geneSymbol, @NotNull String name) {
        if (geneSymbol == null) {
            LOGGER.debug("Skipping extraction for '{}' since gene is missing", name);
            return MutationType.UNKNOWN;
        } else {
            return classifier.determineType(geneSymbol, name);
        }
    }
}
