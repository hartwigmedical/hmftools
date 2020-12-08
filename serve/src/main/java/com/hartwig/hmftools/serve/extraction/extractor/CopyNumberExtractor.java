package com.hartwig.hmftools.serve.extraction.extractor;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CopyNumberExtractor {

    private static final Set<EventType> COPY_NUMBER_EVENTS = Sets.newHashSet(EventType.AMPLIFICATION, EventType.DELETION);

    @NotNull
    private final GeneChecker geneChecker;

    public CopyNumberExtractor(@NotNull final GeneChecker geneChecker) {
        this.geneChecker = geneChecker;
    }

    @Nullable
    public KnownCopyNumber extract(@NotNull String gene, @NotNull EventType type) {
        if (COPY_NUMBER_EVENTS.contains(type) && geneChecker.isValidGene(gene)) {
            return ImmutableKnownCopyNumber.builder().gene(gene).type(toCopyNumberType(type)).build();
        }

        return null;
    }

    @NotNull
    private static CopyNumberType toCopyNumberType(@NotNull EventType eventType) {
        assert COPY_NUMBER_EVENTS.contains(eventType);

        switch (eventType) {
            case AMPLIFICATION:
                return CopyNumberType.AMPLIFICATION;
            case DELETION:
                return CopyNumberType.DELETION;
            default:
                throw new IllegalStateException("Could not convert event type to copy number type: " + eventType);
        }
    }
}
