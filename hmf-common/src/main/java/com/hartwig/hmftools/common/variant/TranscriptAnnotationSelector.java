package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.gene.CanonicalTranscript;
import com.hartwig.hmftools.common.gene.TranscriptRegion;

import org.jetbrains.annotations.NotNull;

public class TranscriptAnnotationSelector {

    @NotNull
    private final Map<String, String> transcripts;

    TranscriptAnnotationSelector(@NotNull final List<CanonicalTranscript> transcripts) {
        this(transcripts.stream().collect(Collectors.toMap(TranscriptRegion::gene, TranscriptRegion::transcriptID)));
    }

    @VisibleForTesting
    TranscriptAnnotationSelector(@NotNull final Map<String, String> geneTranscriptMap) {
        this.transcripts = geneTranscriptMap;
    }

    @NotNull
    public <T extends TranscriptAnnotation> Optional<T> canonical(@NotNull final String gene,
            @NotNull final List<T> annotations) {
        if (transcripts.containsKey(gene)) {
            final String transcript = transcripts.get(gene);
            return annotations.stream().filter(x -> x.transcript().equals(transcript)).findFirst();
        }

        return Optional.empty();
    }
}
