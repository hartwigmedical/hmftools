package com.hartwig.hmftools.sage.read;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

public class ReadContextCounterFactory {

    private final SageConfig config;

    public ReadContextCounterFactory(final SageConfig config) {
        this.config = config;
    }

    public List<ReadContextCounter> create(@NotNull final String sample, @NotNull final List<Candidate> candidates) {
        return candidates.stream()
                .map(x -> new ReadContextCounter(sample, x.variant(), x.readContext(), x.readDepth() < config.maxReadDepthCandidate()))
                .collect(Collectors.toList());
    }

}
