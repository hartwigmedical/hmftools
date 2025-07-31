package com.hartwig.hmftools.sage.append;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.CandidateSerialisation.PRE_v3_5_FLANK_EXTENSION_LENGTH;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.checkGenotypeFields;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.createGenotype;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sage.FragmentLengthCounts;
import com.hartwig.hmftools.sage.quality.BqrRecordMap;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.evidence.FragmentLengthWriter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.phase.AppendVariantPhaser;
import com.hartwig.hmftools.sage.pipeline.EvidenceStage;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class RegionAppendTask implements Callable<Void>
{
    private final ChrBaseRegion mRegion;
    private final int mTaskId;

    private final SageAppendConfig mConfig;
    private final EvidenceStage mEvidenceStage;
    private final IndexedFastaSequenceFile mRefGenomeFile;
    private final RefGenomeSource mRefGenome;
    private final FragmentLengthWriter mFragmentLengths;
    private final AppendVariantPhaser mVariantPhaser;
    private final SamSlicerFactory mSamSlicerFactory;

    private final List<VariantContext> mOriginalVariants;
    private final List<VariantContext> mFinalVariants;

    public RegionAppendTask(
            final int taskId, final ChrBaseRegion region, final List<VariantContext> variants,
            final SageAppendConfig config, final IndexedFastaSequenceFile refGenome, final Map<String, BqrRecordMap> qualityRecalibrationMap,
            final FragmentLengthWriter fragmentLengths, final MsiJitterCalcs msiJitterCalcs)
    {
        mTaskId = taskId;
        mRegion = region;
        mOriginalVariants = variants;
        mFragmentLengths = fragmentLengths;

        mFinalVariants = Lists.newArrayList();

        mConfig = config;
        mRefGenomeFile = refGenome;
        mRefGenome = new RefGenomeSource(mRefGenomeFile);

        mSamSlicerFactory = new SamSlicerFactory();
        mSamSlicerFactory.buildBamReaders(Collections.emptyList(), Collections.emptyList(), mConfig.Common, mRefGenomeFile);

        mVariantPhaser = new AppendVariantPhaser();

        mEvidenceStage = new EvidenceStage(
                config.Common, mRefGenome, qualityRecalibrationMap, msiJitterCalcs, mVariantPhaser, mSamSlicerFactory);
    }

    public List<VariantContext> finalVariants() { return mFinalVariants; }

    @Override
    public Void call()
    {
        SG_LOGGER.trace("{}: region({}) finding evidence", mTaskId, mRegion);

        ChrBaseRegion extendedRegion = new ChrBaseRegion(
                mRegion.Chromosome,
                max(1, mRegion.start() - PRE_v3_5_FLANK_EXTENSION_LENGTH),
                mRegion.end() + PRE_v3_5_FLANK_EXTENSION_LENGTH);

        RefSequence refSequence = new RefSequence(extendedRegion, mRefGenome);

        List<Candidate> candidates = mOriginalVariants.stream()
                .map(x -> CandidateSerialisation.toCandidate(x, refSequence))
                .filter(x -> x != null)
                .collect(Collectors.toList());

        if(candidates.size() < mOriginalVariants.size())
        {
            SG_LOGGER.error("region({}) failed to recreate variant context for {} variants",
                    mRegion, mOriginalVariants.size() - candidates.size());
            System.exit(1);
        }

        mVariantPhaser.registerLocalPhaseSets(candidates, mOriginalVariants);

        ReadContextCounters readContextCounters = mEvidenceStage.findEvidence
                (mRegion, "reference", mConfig.Common.ReferenceIds, candidates, mConfig.Common.ReferenceIds);

        createFinalVariants(readContextCounters, mConfig.Common.ReferenceIds);

        mVariantPhaser.populateLocalPhaseSetInfo(candidates, mFinalVariants);

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
                    FragmentLengthCounts fragmentLengthData = sampleCounters.get(s).fragmentLengthCounts();
                    mFragmentLengths.writeVariantFragmentLength(variantInfo, sampleId, fragmentLengthData);
                }
            }
        }

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        mSamSlicerFactory.closeSamReaders();

        return null;
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
            final VariantContext variantContext, final List<ReadContextCounter> readCounters, final List<String> sampleIds)
    {
        VariantContextBuilder builder = new VariantContextBuilder(variantContext);

        List<Genotype> genotypes = Lists.newArrayList();

        for(Genotype genotype : variantContext.getGenotypes())
        {
            genotypes.add(checkGenotypeFields(genotype));
        }

        for(int i = 0; i < readCounters.size(); ++i)
        {
            ReadContextCounter readCounter = readCounters.get(i);
            Genotype genotype = createGenotype(readCounter, sampleIds.get(i));
            genotypes.add(genotype);
        }

        return builder.genotypes(genotypes).make();
    }
}
