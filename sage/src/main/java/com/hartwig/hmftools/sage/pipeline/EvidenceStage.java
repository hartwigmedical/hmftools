package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextEvidence;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class EvidenceStage
{
    private final ReadContextEvidence mReadContextEvidence;

    public EvidenceStage(
            final SageConfig config, final ReferenceSequenceFile refGenome,
            final Map<String, QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mReadContextEvidence = new ReadContextEvidence(config, refGenome, qualityRecalibrationMap);
    }

    public CompletableFuture<ReadContextCounters> findEvidence(
            final ChrBaseRegion region, final String sampleType, final List<String> samples, final List<String> sampleBams,
            final CompletableFuture<List<Candidate>> candidates)
    {
        // Scan tumors for evidence
        return candidates.thenCompose(initialCandidates ->
        {
            final String primarySample = samples.isEmpty() ? "PRIMARY" : samples.get(0);

            final ReadContextCounters result = new ReadContextCounters(primarySample, initialCandidates);

            // SG_LOGGER.trace("region({}) gathering {} evidence for {} candidates", region, sampleType, initialCandidates.size());

            CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
            for(int i = 0; i < samples.size(); i++)
            {
                final String sample = samples.get(i);
                final String sampleBam = sampleBams.get(i);

                // SG_LOGGER.trace("region({}) tumor sample({}) gathering evidence for {} candidates", region, sample, initialCandidates.size());

                done = done.<List<ReadContextCounter>>thenApply(x -> mReadContextEvidence.get(initialCandidates, sample, sampleBam)).thenAccept(result::addCounters);
            }

            SG_LOGGER.trace("region({}) gathered {} evidence for {} variants", region, sampleType, result.readContextCounters().size());

            return done.thenApply(x -> result);
        });
    }
}
