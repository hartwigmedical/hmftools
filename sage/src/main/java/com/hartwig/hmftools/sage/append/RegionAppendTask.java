package com.hartwig.hmftools.sage.append;

import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.createGenotype;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.evidence.FragmentLengthData;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.EvidenceStage;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class RegionAppendTask implements Callable
{
    private final ChrBaseRegion mRegion;
    private final int mTaskId;

    private final SageAppendConfig mConfig;
    private final EvidenceStage mEvidenceStage;
    private final IndexedFastaSequenceFile mRefGenomeFile;
    private final RefGenomeSource mRefGenome;
    private final FragmentLengths mFragmentLengths;

    private final List<VariantContext> mOriginalVariants;
    private final List<VariantContext> mFinalVariants;

    public RegionAppendTask(
            final int taskId, final ChrBaseRegion region, final List<VariantContext> variants,
            final SageAppendConfig config, final IndexedFastaSequenceFile refGenome,
            final Map<String, BqrRecordMap> qualityRecalibrationMap, final FragmentLengths fragmentLengths)
    {
        mTaskId = taskId;
        mRegion = region;
        mOriginalVariants = variants;
        mFragmentLengths = fragmentLengths;

        mFinalVariants = Lists.newArrayList();

        mConfig = config;
        mRefGenomeFile = refGenome;
        mRefGenome = new RefGenomeSource(mRefGenomeFile);

        SamSlicerFactory samSlicerFactory = new SamSlicerFactory();
        samSlicerFactory.buildBamReaders(Collections.emptyList(), Collections.emptyList(), mConfig.Common, mRefGenomeFile);

        mEvidenceStage = new EvidenceStage(config.Common, mRefGenome, qualityRecalibrationMap, new PhaseSetCounter(), samSlicerFactory);
    }

    public List<VariantContext> finalVariants() { return mFinalVariants; }

    @Override
    public Long call()
    {
        SG_LOGGER.trace("{}: region({}) finding evidence", mTaskId, mRegion);

        RefSequence refSequence = new RefSequence(mRegion, mRefGenomeFile);

        List<Candidate> candidates = mOriginalVariants.stream()
                .map(x -> CandidateSerialization.toCandidate(x, refSequence)).collect(Collectors.toList());

        ReadContextCounters readContextCounters = mEvidenceStage.findEvidence
                (mRegion, "reference", mConfig.Common.ReferenceIds, candidates, false);

        createFinalVariants(readContextCounters, mConfig.Common.ReferenceIds);

        if(mConfig.Common.WriteFragmentLengths)
        {
            for(int i = 0; i < mOriginalVariants.size(); ++i)
            {
                Candidate variant = candidates.get(i);

                String variantInfo = format("%s:%d %s>%s",
                        variant.chromosome(), variant.position(), variant.variant().ref(), variant.variant().alt());

                List<ReadContextCounter> sampleCounters = readContextCounters.getReadCounters(i);

                for(int s = 0; s < mConfig.Common.ReferenceIds.size(); ++s)
                {
                    String sampleId = mConfig.Common.ReferenceIds.get(s);
                    FragmentLengthData fragmentLengthData = sampleCounters.get(s).fragmentLengths();
                    mFragmentLengths.writeVariantFragmentLength(variantInfo, sampleId, fragmentLengthData);
                }
            }
        }

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        return (long)0;
    }

    public void createFinalVariants(final ReadContextCounters readContextCounters, final List<String> sampleIds)
    {
        for(int i = 0; i < mOriginalVariants.size(); ++i)
        {
            VariantContext origVariant = mOriginalVariants.get(i);

            List<ReadContextCounter> sampleCounters = readContextCounters.getReadCounters(i);
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
