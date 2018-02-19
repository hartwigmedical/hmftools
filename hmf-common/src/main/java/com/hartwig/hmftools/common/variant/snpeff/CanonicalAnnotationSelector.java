package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.gene.CanonicalTranscript;
import com.hartwig.hmftools.common.gene.TranscriptRegion;

import org.jetbrains.annotations.NotNull;

public class CanonicalAnnotationSelector {

    private final Map<String, String> transcripts;

    public CanonicalAnnotationSelector(final List<CanonicalTranscript> transcripts) {
        this(transcripts.stream().collect(Collectors.toMap(TranscriptRegion::gene, TranscriptRegion::transcriptID)));
    }

    @VisibleForTesting
    CanonicalAnnotationSelector(final Map<String, String> geneTranscriptMap) {
        this.transcripts = geneTranscriptMap;
    }

    @NotNull
    public Optional<VariantAnnotation> canonical(@NotNull final String gene, @NotNull final List<VariantAnnotation> annotations) {
        if (transcripts.containsKey(gene)) {
            final String transcript = transcripts.get(gene);
            return annotations.stream().filter(x -> x.featureID().equals(transcript)).findFirst();
        }

        return Optional.empty();
    }

}
