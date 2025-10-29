package com.hartwig.hmftools.sage.append;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.CONSENSUS_TAG_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MIN_COORDS_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.CandidateSerialisation.PRE_v3_5_FLANK_EXTENSION_LENGTH;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.createGenotype;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_SEQ_TECH_BASE_QUAL;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
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
import htsjdk.variant.variantcontext.GenotypeBuilder;
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

        // first gather existing genotype attributes from existing genotypes
        Set<String> existingGenotypeTags = Sets.newHashSet();

        for(Genotype genotype : variantContext.getGenotypes())
        {
            genotype.getExtendedAttributes().keySet().forEach(x -> existingGenotypeTags.add(x));
        }

        Set<String> expectedGenotypeTags = Sets.newHashSet();

        List<Genotype> newGenotypes = Lists.newArrayList();

        for(int i = 0; i < readCounters.size(); ++i)
        {
            ReadContextCounter readCounter = readCounters.get(i);
            Genotype genotype = createGenotype(readCounter, sampleIds.get(i), existingGenotypeTags);
            newGenotypes.add(genotype);

            if(expectedGenotypeTags.isEmpty())
                expectedGenotypeTags.addAll(genotype.getExtendedAttributes().keySet());
        }

        List<Genotype> genotypes = Lists.newArrayList();

        // check that existing genotypes have an entry for all tags
        for(Genotype genotype : variantContext.getGenotypes())
        {
            genotypes.add(checkGenotypeFields(genotype, expectedGenotypeTags));
        }

        newGenotypes.forEach(x -> genotypes.add(x));

        return builder.genotypes(genotypes).make();
    }

    // built in fields: GT:AD:DP, not populated: RAD:RDP
    // v4.0 attributes: ABQ:    AF:AMBQ:AMMQ:AMQ:         RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RSB:SAC:SB:UMI_CNT
    // v4.1 attributes: ABQ:AED:AF:AMBQ:AMMQ:AMQ:MUC:RABQ:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RSB:SAC:SB:UMI_CNT

    // added in v4.1: MIN_COORDS_FLAG, AVG_RAW_BASE_QUAL, AVG_READ_EDGE_DISTANCE
    // added in v4.2: AVG_EDGE_DISTANCE_PERC, replacing AVG_READ_EDGE_DISTANCE

    private static Genotype checkGenotypeFields(final Genotype genotype, final Set<String> expectedTags)
    {
        // checks that existing genotype contain all required fields (eg append is adding any new ones) and if not adds defaults
        if(expectedTags.stream().allMatch(x -> genotype.getExtendedAttributes().containsKey(x)))
            return genotype;

        // otherwise add meaningful defaults where possible
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);

        for(String vcfTag : expectedTags)
        {
            if(genotype.getExtendedAttributes().containsKey(vcfTag))
                continue;

            if(vcfTag.equals(MIN_COORDS_COUNT))
                genotypeBuilder.attribute(MIN_COORDS_COUNT, 0);
            else if(vcfTag.equals(AVG_SEQ_TECH_BASE_QUAL))
                genotypeBuilder.attribute(AVG_SEQ_TECH_BASE_QUAL, new int[] {0, 0});
            else if(vcfTag.equals(AVG_EDGE_DISTANCE_PERC))
                genotypeBuilder.attribute(AVG_EDGE_DISTANCE_PERC, new double[] {0, 0});
            else if(vcfTag.equals(UMI_TYPE_COUNTS))
                genotypeBuilder.attribute(UMI_TYPE_COUNTS, new int[CONSENSUS_TAG_TYPE_COUNT]);
            else
                genotypeBuilder.attribute(vcfTag, "");
        }

        return genotypeBuilder.make();
    }
}
