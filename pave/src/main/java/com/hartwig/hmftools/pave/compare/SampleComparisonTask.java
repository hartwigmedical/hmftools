package com.hartwig.hmftools.pave.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.compare.RefVariantData.loadVariantsFromVcf;
import static com.hartwig.hmftools.pave.impact.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.impact.PaveUtils.findVariantImpacts;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.hasCodingEffectDiff;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.CODING_EFFECT;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.HGVS_CODING;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.HGVS_PROTEIN;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.REPORTED;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.annotation.Reportability;
import com.hartwig.hmftools.pave.VariantData;
import com.hartwig.hmftools.pave.impact.VariantImpactBuilder;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

public class SampleComparisonTask implements Callable
{
    private final int mTaskId;
    private final ComparisonConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;
    private final GeneDataCache mGeneDataCache;
    private final Reportability mReportability;

    private final ComparisonWriter mWriter;

    private final List<String> mSampleIds;
    private final Map<String,List<RefVariantData>> mSampleVariantsCache;

    private static final String PC_PROCESS = "Process";
    private final Map<String,PerformanceCounter> mPerfCounters;

    private int mTotalComparisons;
    private int mMatchedCount;

    public SampleComparisonTask(
            int taskId, final ComparisonConfig config, final RefGenomeInterface refGenome,
            final ComparisonWriter writer, final GeneDataCache geneDataCache, final Map<String,List<RefVariantData>> sampleVariantsCache)
    {
        mTaskId = taskId;
        mGeneDataCache = geneDataCache;
        mConfig = config;
        mWriter = writer;
        mSampleVariantsCache = sampleVariantsCache;

        mSampleIds = Lists.newArrayList();
        mImpactClassifier = new ImpactClassifier(refGenome);
        mImpactBuilder = new VariantImpactBuilder(mGeneDataCache);
        mReportability = new Reportability(mGeneDataCache.getDriverPanel());

        mTotalComparisons = 0;
        mMatchedCount = 0;

        mPerfCounters = Maps.newHashMap();
        mPerfCounters.put(PC_PROCESS, new PerformanceCounter(PC_PROCESS));
    }

    public List<String> getSampleIds() { return mSampleIds; }
    public int totalComparisons() { return mTotalComparisons; }
    public int matchedCount() { return mMatchedCount; }

    public Map<String,PerformanceCounter> getPerfCounters() { return mPerfCounters; }

    @Override
    public Long call()
    {
        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            String sampleId = mSampleIds.get(i);

            checkSampleDiffs(sampleId);

            if(i > 0 && (i % 100) == 0)
            {
                PV_LOGGER.info("{}: processed {} samples", mTaskId, i);
            }
        }

        PV_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());

        return (long)0;
    }

    private void checkSampleDiffs(final String sampleId)
    {
        List<RefVariantData> refVariants = null;

        if(!mSampleVariantsCache.isEmpty())
        {
            refVariants = mSampleVariantsCache.get(sampleId);
        }
        else
        {
            String sampleVcf = convertWildcardSamplePath(mConfig.SampleVCF, sampleId);
            refVariants = loadVariantsFromVcf(sampleVcf, mConfig);
        }

        if(refVariants == null || refVariants.isEmpty())
        {
            PV_LOGGER.error("failed to load sample({}) variants", sampleId);
            return;
        }

        mPerfCounters.get(PC_PROCESS).start();

        for(RefVariantData refVariant :refVariants)
        {
            // generate variant impact data and then write comparison results to TSV

            boolean isDriverGene = mGeneDataCache.isDriverPanelGene(refVariant.Gene);

            if(mConfig.OnlyDriverGenes && !isDriverGene)
                continue;

            VariantData variant = new VariantData(
                    refVariant.Chromosome, refVariant.Position, refVariant.Ref, refVariant.Alt);

            variant.setVariantDetails(refVariant.LocalPhaseSet, refVariant.Microhomology, refVariant.RepeatSequence, refVariant.RepeatCount);

            try
            {
                variant.setRealignedVariant(createRightAlignedVariant(variant, mImpactClassifier.refGenome()));

                findVariantImpacts(variant, mImpactClassifier, mGeneDataCache);

                processPhasedVariants(variant.localPhaseSet(), sampleId, refVariants);

                if(!variant.hasLocalPhaseSet())
                    processVariant(sampleId, variant, refVariant, isDriverGene);
            }
            catch(Exception e)
            {
                PV_LOGGER.error("error processing var({}): {}", variant, e.toString());
                e.printStackTrace();
            }
        }

        processPhasedVariants(NO_LOCAL_PHASE_SET, sampleId, refVariants);
        mImpactClassifier.phasedVariants().clear();

        mPerfCounters.get(PC_PROCESS).stop();

        PV_LOGGER.debug("{}: sample({}) processed {} variants", mTaskId, sampleId, refVariants.size());
    }

    private void processPhasedVariants(int currentLocalPhaseSet, final String sampleId, final List<RefVariantData> refVariants)
    {
        List<VariantData> variants = mImpactClassifier.processPhasedVariants(currentLocalPhaseSet);

        if(variants != null)
        {
            for(VariantData variant : variants)
            {
                RefVariantData refVariant = refVariants.stream().filter(x -> x.matches(variant)).findFirst().orElse(null);
                boolean isDriverGene = mGeneDataCache.isDriverPanelGene(refVariant.Gene);
                processVariant(sampleId, variant, refVariant, isDriverGene);
            }
        }
    }

    private void processVariant(
            final String sampleId, final VariantData variant, final RefVariantData refVariant, boolean isDriverGene)
    {
        ++mTotalComparisons;

        if(refVariant.Gene.isEmpty() && variant.getImpacts().isEmpty())
        {
            ++mMatchedCount;
            return;
        }

        VariantImpact variantImpact = mImpactBuilder.createVariantImpact(variant);

        boolean reportable = isReported(variant, variantImpact, refVariant);

        if(reportable)
            variant.markReported();

        List<String> diffs = Lists.newArrayList();

        if(variantImpact == null)
        {
            // hasDiff = refVariant.Reported;
            return;
        }
        else
        {
            if(mConfig.checkDiffType(REPORTED) && refVariant.Reported != variant.reported())
                diffs.add(format("Reported(%s/%s)", variant.reported(), refVariant.Reported));

            if(mConfig.checkDiffType(CODING_EFFECT))
            {
                if(!variantImpact.CanonicalEffect.equals(refVariant.CanonicalEffect))
                    diffs.add(format("CanonicalEffect(%s/%s)", variantImpact.CanonicalEffect, refVariant.CanonicalEffect));

                if(!variantImpact.WorstCodingEffect.equals(refVariant.WorstCodingEffect))
                    diffs.add(format("WorstCodingEffect(%s/%s)", variantImpact.WorstCodingEffect, refVariant.WorstCodingEffect));

                if(hasCodingEffectDiff(variantImpact, refVariant))
                    diffs.add(format("CodingEffect(%s/%s)", variantImpact.CanonicalCodingEffect, refVariant.CanonicalCodingEffect));
            }

            VariantTransImpact transImpact = variant.getCanonicalTransImpacts(refVariant.Gene);

            if(transImpact != null)
            {
                VariantTransImpact raTransImpact = variant.getRealignedImpact(refVariant.Gene, transImpact);

                if(mConfig.checkDiffType(HGVS_CODING))
                {
                    //if(hasHgvsCodingDiff(transImpact, raTransImpact, refVariant))
                    if(!variantImpact.CanonicalHgvsCoding.equals(refVariant.HgvsCodingImpact))
                        diffs.add(format("HgvsCoding(%s/%s)", variantImpact.CanonicalHgvsCoding, refVariant.HgvsCodingImpact));
                }

                if(mConfig.checkDiffType(HGVS_PROTEIN))
                {
                    //if(hasHgvsProteinDiff(variant, transImpact, raTransImpact, refVariant))
                    if(!variantImpact.CanonicalHgvsProtein.equals(refVariant.HgvsProteinImpact))
                        diffs.add(format("HgvsProtein(%s/%s)", variantImpact.CanonicalHgvsProtein, refVariant.HgvsProteinImpact));
                }
            }
        }

        if(diffs.isEmpty())
        {
            ++mMatchedCount;
        }
        else
        {
            logComparison(sampleId, refVariant, variant, variantImpact);
        }

        if(!diffs.isEmpty() || mConfig.WriteMatches)
            mWriter.writeVariantDiff(sampleId, variant, variantImpact, refVariant, isDriverGene, diffs);
    }

    private boolean isReported(final VariantData variant, final VariantImpact variantImpact, final RefVariantData refVariant)
    {
        return mReportability.isReported(variant, variantImpact, refVariant.IsHotspot);
    }

    private void logComparison(
            final String sampleId, final RefVariantData refVariant, final VariantData variant, final VariantImpact variantImpact)
    {
        if(!PV_LOGGER.isDebugEnabled())
            return;

        if(variant.getImpacts().isEmpty())
        {
            PV_LOGGER.trace("sample({}) var({}) no gene impacts found vs snpEff({} : {})",
                    sampleId, variant.toString(), refVariant.Gene, refVariant.CanonicalCodingEffect);
            return;
        }

        if(!variant.getImpacts().containsKey(refVariant.Gene))
        {
            PV_LOGGER.trace("sample({}) var({}) diff gene and canonical coding: pave({} : {}) snpEff({} : {})",
                    sampleId, variant.toString(), variantImpact.GeneName, variantImpact.CanonicalCodingEffect,
                    refVariant.Gene, refVariant.CanonicalCodingEffect);
            return;
        }

        PV_LOGGER.trace("sample({}) var({}) diff canonical coding: pave({} : {}) snpEff({} : {})",
                sampleId, variant.toString(), variantImpact.GeneName, variantImpact.CanonicalCodingEffect,
                refVariant.Gene, refVariant.CanonicalCodingEffect);
    }
}
