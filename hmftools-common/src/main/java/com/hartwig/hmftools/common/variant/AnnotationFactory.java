package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class AnnotationFactory {

    private static final Logger LOGGER = LogManager.getLogger(AnnotationFactory.class);

    private static final String START_IDENTIFIER = "ANN=";
    private static final String END_IDENTIFIER = ";";
    private static final String ANNOTATION_SEPARATOR = ",";
    private static final String FIELD_SEPARATOR = "\\|";

    private static final int EXPECTED_FIELD_SIZE_PER_ANNOTATION = 16;

    private AnnotationFactory() {
    }

    @NotNull
    static List<Annotation> fromVCFInfoField(@NotNull final String info) {
        final List<Annotation> annotations = Lists.newArrayList();
        final int startIndex = info.indexOf(START_IDENTIFIER);
        if (startIndex >= 0) {
            String fullAnnotationString = info.substring(startIndex + START_IDENTIFIER.length());
            final int endIndex = fullAnnotationString.indexOf(END_IDENTIFIER);
            if (endIndex > 0) {
                fullAnnotationString = fullAnnotationString.substring(0, endIndex - 1);
            }
            for (final String annotationString : fullAnnotationString.split(ANNOTATION_SEPARATOR)) {
                final String[] parts = enforceMinLength(annotationString.split(FIELD_SEPARATOR),
                        EXPECTED_FIELD_SIZE_PER_ANNOTATION);
                if (parts.length == EXPECTED_FIELD_SIZE_PER_ANNOTATION) {
                    annotations.add(
                            new Annotation(parts[0], toConsequence(parts[1]), parts[2], parts[3], parts[4], parts[5],
                                    parts[6], parts[7], parts[8], parts[9], parts[10], parts[11], parts[12], parts[13],
                                    parts[14], parts[15]));
                } else {
                    LOGGER.warn("Annotation found with invalid field count: " + annotationString);
                }
            }

        }
        return annotations;
    }

    @NotNull
    private static String[] enforceMinLength(@NotNull String[] parts, int minSize) {
        if (parts.length > minSize) {
            return parts;
        } else {
            final String[] values = new String[minSize];
            for (int i = 0; i < minSize; i++) {
                values[i] = i < parts.length ? parts[i] : Strings.EMPTY;
            }
            System.arraycopy(parts, 0, values, 0, parts.length);

            return values;
        }
    }

    @NotNull
    private static VariantConsequence toConsequence(@NotNull final String annotation) {
        for (final VariantConsequence consequence : VariantConsequence.values()) {
            if (consequence.sequenceOntologyTerm().equals(annotation)) {
                return consequence;
            }
        }
        return VariantConsequence.OTHER;
    }
}
