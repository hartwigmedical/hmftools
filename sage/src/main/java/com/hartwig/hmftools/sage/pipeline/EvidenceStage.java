package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;

import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextEvidence;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContextCounters;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class EvidenceStage
{
    private final ReadContextEvidence mReadContextEvidence;

    public EvidenceStage(@NotNull final SageConfig config, @NotNull final ReferenceSequenceFile refGenome,
            @NotNull final Map<String, QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mReadContextEvidence = new ReadContextEvidence(config, refGenome, qualityRecalibrationMap);
    }

    @NotNull
    public CompletableFuture<ReadContextCounters> evidence(@NotNull final List<String> samples, @NotNull final List<String> sampleBams,
            @NotNull final CompletableFuture<List<Candidate>> candidates)
    {
        // Scan tumors for evidence
        return candidates.thenCompose(initialCandidates ->
        {
            final String primarySample = samples.isEmpty() ? "PRIMARY" : samples.get(0);

            final ReadContextCounters result = new ReadContextCounters(primarySample, initialCandidates);

            CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
            for(int i = 0; i < samples.size(); i++)
            {
                final String sample = samples.get(i);
                final String sampleBam = sampleBams.get(i);

                done = done.thenApply(x -> mReadContextEvidence.get(initialCandidates, sample, sampleBam)).thenAccept(result::addCounters);
            }

            return done.thenApply(x -> result);
        });
    }

}
