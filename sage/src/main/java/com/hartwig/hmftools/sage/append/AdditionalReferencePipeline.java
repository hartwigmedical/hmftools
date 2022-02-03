package com.hartwig.hmftools.sage.append;

import static java.util.concurrent.CompletableFuture.supplyAsync;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.createGenotype;

import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.EvidenceStage;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class AdditionalReferencePipeline
{
    private final SageConfig mConfig;
    private final EvidenceStage mEvidenceStage;
    private final ReferenceSequenceFile mRefGenome;
    private final Executor mExecutor;

    public AdditionalReferencePipeline(
            final SageConfig config, final Executor executor, ReferenceSequenceFile refGenome,
            final Map<String, QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mEvidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap, new PhaseSetCounter());
        mExecutor = executor;
    }

    public CompletableFuture<List<VariantContext>> appendReference(final ChrBaseRegion region, final List<VariantContext> variants)
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

        // final CompletableFuture<ReadContextCounters> evidenceFutures =
        //        mEvidenceStage.findEvidenceOld(region, "reference", mConfig.ReferenceIds, mConfig.ReferenceBams, candidateFutures, false);

        return null; // evidenceFutures.thenApply(x -> update(x, variants));
    }

    public List<VariantContext> update(
            final ReadContextCounters readContextCounters, final List<VariantContext> variantContexts, final List<String> sampleIds)
    {
        final List<VariantContext> result = Lists.newArrayList();
        for(VariantContext old : variantContexts)
        {
            VariantHotspot variant = CandidateSerialization.toVariantHotspot(old);
            List<ReadContextCounter> counters = readContextCounters.getVariantReadCounters(variant);
            result.add(addGenotype(old, counters, sampleIds));
        }

        return result;
    }

    private static VariantContext addGenotype(
            final VariantContext parent, final List<ReadContextCounter> readCounters, final List<String> sampleIds)
    {
        final VariantContextBuilder builder = new VariantContextBuilder(parent);
        final List<Genotype> genotypes = Lists.newArrayList(parent.getGenotypes());

        for(int i = 0; i < readCounters.size(); ++i)
        {
            ReadContextCounter readCounter = readCounters.get(i);
            Genotype genotype = createGenotype(readCounter, sampleIds.get(i));
            genotypes.add(genotype);
        }
        return builder.genotypes(genotypes).make();
    }

}
