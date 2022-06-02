package com.hartwig.hmftools.sage.append;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.createGenotype;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.EvidenceStage;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class RegionAppendTask implements Callable
{
    private final ChrBaseRegion mRegion;
    private final int mTaskId;

    private final SageConfig mConfig;
    private final EvidenceStage mEvidenceStage;
    private final IndexedFastaSequenceFile mRefGenomeFile;
    private final RefGenomeSource mRefGenome;

    private final List<VariantContext> mOriginalVariants;
    private final List<VariantContext> mFinalVariants;

    public RegionAppendTask(
            final int taskId, final ChrBaseRegion region, final List<VariantContext> variants,
            final SageConfig config, final IndexedFastaSequenceFile refGenome,
            final Map<String,QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mTaskId = taskId;
        mRegion = region;
        mOriginalVariants = variants;
        mFinalVariants = Lists.newArrayList();

        mConfig = config;
        mRefGenomeFile = refGenome;
        mRefGenome = new RefGenomeSource(mRefGenomeFile);

        SamSlicerFactory samSlicerFactory = new SamSlicerFactory();
        samSlicerFactory.buildBamReaders(mConfig, mRefGenomeFile);

        mEvidenceStage = new EvidenceStage(config, mRefGenome, qualityRecalibrationMap, new PhaseSetCounter(), samSlicerFactory);
    }

    public List<VariantContext> finalVariants() { return mFinalVariants; }

    @Override
    public Long call()
    {
        SG_LOGGER.trace("{}: region({}) finding evidence", mTaskId, mRegion);

        RefSequence refSequence = new RefSequence(mRegion, mRefGenomeFile);

        List<Candidate> candidates = mOriginalVariants.stream()
                .map(x -> CandidateSerialization.toCandidate(x, refSequence)).collect(Collectors.toList());

        ReadContextCounters normalEvidence = mEvidenceStage.findEvidence
                (mRegion, "reference", mConfig.ReferenceIds, candidates, false);

        createFinalVariants(normalEvidence, mConfig.ReferenceIds);

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        return (long)0;
    }

    public void createFinalVariants(final ReadContextCounters readContextCounters, final List<String> sampleIds)
    {
        for(int i = 0; i < mOriginalVariants.size(); ++i)
        {
            VariantContext origVariant = mOriginalVariants.get(i);

            List<ReadContextCounter> sampleCounters = readContextCounters.getReadCounters(i); // getVariantReadCounters(variant);
            mFinalVariants.add(addGenotype(origVariant, sampleCounters, sampleIds));
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
