package com.hartwig.hmftools.sage.append;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.createGenotype;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.EvidenceStage;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class RegionAppendTask implements Callable
{
    private final ChrBaseRegion mRegion;
    private final int mTaskId;

    private final SageConfig mConfig;
    private final EvidenceStage mEvidenceStage;
    private final ReferenceSequenceFile mRefGenome;

    private final List<VariantContext> mOriginalVariants;
    private final List<VariantContext> mFinalVariants;

    public RegionAppendTask(
            final int taskId, final ChrBaseRegion region, final List<VariantContext> variants,
            final SageConfig config, final ReferenceSequenceFile refGenome,
            final Map<String,QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mTaskId = taskId;
        mRegion = region;
        mOriginalVariants = variants;
        mFinalVariants = Lists.newArrayList();

        mConfig = config;
        mRefGenome = refGenome;
        mEvidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap, new PhaseSetCounter());
    }

    public List<VariantContext> finalVariants() { return mFinalVariants; }

    @Override
    public Long call()
    {
        SG_LOGGER.trace("{}: region({}) finding evidence", mTaskId, mRegion);

        RefSequence refSequence = new RefSequence(mRegion, mRefGenome);

        List<Candidate> candidates = mOriginalVariants.stream()
                .map(x -> CandidateSerialization.toCandidate(x, refSequence)).collect(Collectors.toList());

        ReadContextCounters normalEvidence = mEvidenceStage.findEvidence
                (mRegion, "reference", mConfig.ReferenceIds, mConfig.ReferenceBams, candidates, false);

        createFinalVariants(normalEvidence, mConfig.ReferenceIds);

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        return (long)0;
    }

    public void createFinalVariants(final ReadContextCounters readContextCounters, final List<String> sampleIds)
    {
        for(VariantContext origVariant : mOriginalVariants)
        {
            VariantHotspot variant = CandidateSerialization.toVariantHotspot(origVariant);
            List<ReadContextCounter> counters = readContextCounters.getVariantReadCounters(variant);
            mFinalVariants.add(addGenotype(origVariant, counters, sampleIds));
        }
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
