package com.hartwig.hmftools.sage.pipeline;

import static java.util.concurrent.CompletableFuture.supplyAsync;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.CandidateSerialization;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextCounters;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class AdditionalReferencePipeline
{
    private final SageConfig mConfig;
    private final EvidenceStage mEvidenceStage;
    private final ReferenceSequenceFile mRefGenome;
    private final Executor mExecutor;

    public AdditionalReferencePipeline(@NotNull final SageConfig config, @NotNull final Executor executor, ReferenceSequenceFile refGenome,
            @NotNull final Map<String, QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mEvidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap);
        mExecutor = executor;
    }

    @NotNull
    public CompletableFuture<List<VariantContext>> appendReference(final BaseRegion region, final List<VariantContext> variants)
    {
        if(variants.isEmpty())
        {
            return CompletableFuture.completedFuture(Lists.newArrayList());
        }

        final CompletableFuture<RefSequence> refSequenceFuture = supplyAsync(() ->
        {
            if(region.start() == 1)
            {
                SG_LOGGER.info("Processing chromosome {}", region.chromosome());
            }
            return new RefSequence(region, mRefGenome);
        }, mExecutor);

        final CompletableFuture<List<Candidate>> candidateFutures = refSequenceFuture.thenApply(x -> variants.stream()
                .map(y -> CandidateSerialization.toCandidate(y, x))
                .collect(Collectors.toList()));

        final CompletableFuture<ReadContextCounters> evidenceFutures =
                mEvidenceStage.findEvidence(mConfig.ReferenceIds, mConfig.ReferenceBams, candidateFutures);

        return evidenceFutures.thenApply(x -> update(x, variants));
    }

    public List<VariantContext> update(final ReadContextCounters readContextCounters, final List<VariantContext> variantContexts)
    {
        final List<VariantContext> result = Lists.newArrayList();
        for(VariantContext old : variantContexts)
        {
            VariantHotspot variant = CandidateSerialization.toVariantHotspot(old);
            List<ReadContextCounter> counters = readContextCounters.readContextCounters(variant);
            result.add(SageVariantContextFactory.addGenotype(old, counters));
        }

        return result;
    }
}
