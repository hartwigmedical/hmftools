package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.pave.PaveApplication.findVariantImpacts;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.hasCodingEffectDiff;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.hasHgvsCodingDiff;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.hasHgvsProteinDiff;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.CODING_EFFECT;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.HGVS_CODING;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.HGVS_PROTEIN;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.REPORTED;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.ImpactClassifier;
import com.hartwig.hmftools.pave.VariantData;
import com.hartwig.hmftools.pave.VariantImpactBuilder;
import com.hartwig.hmftools.pave.VariantTransImpact;

public class SampleComparisonTask implements Callable
{
    private final int mTaskId;
    private final ComparisonConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;
    private final GeneDataCache mGeneDataCache;
    private final ReportablePredicate mReportableOncoGenes;
    private final ReportablePredicate mReportableTsgGenes;

    private final DatabaseAccess mDbAccess;
    private final ComparisonWriter mWriter;

    private final List<String> mSampleIds;
    private final Map<String,List<RefVariantData>> mSampleVariantsCache;
    private final PerformanceCounter mPerfCounter;

    private int mTotalComparisons;
    private int mMatchedCount;

    public SampleComparisonTask(
            int taskId, final ComparisonConfig config, RefGenomeInterface refGenome, final DatabaseAccess dbAccess,
            final ComparisonWriter writer, final GeneDataCache geneDataCache, final Map<String,List<RefVariantData>> sampleVariantsCache)
    {
        mTaskId = taskId;
        mGeneDataCache = geneDataCache;
        mConfig = config;
        mDbAccess = dbAccess;
        mWriter = writer;
        mSampleVariantsCache = sampleVariantsCache;

        mSampleIds = Lists.newArrayList();
        mImpactClassifier = new ImpactClassifier(refGenome);
        mImpactBuilder = new VariantImpactBuilder(mGeneDataCache);
        mReportableOncoGenes = new ReportablePredicate(DriverCategory.ONCO, mGeneDataCache.getDriverPanel());
        mReportableTsgGenes = new ReportablePredicate(DriverCategory.TSG, mGeneDataCache.getDriverPanel());

        mPerfCounter = new PerformanceCounter("ClassifyVariant");
        mTotalComparisons = 0;
        mMatchedCount = 0;
    }

    public List<String> getSampleIds() { return mSampleIds; }
    public int totalComparisons() { return mTotalComparisons; }
    public int matchedCount() { return mMatchedCount; }

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
        List<RefVariantData> refVariants = mDbAccess != null ?
                DataLoader.loadSampleDatabaseRecords(sampleId, mDbAccess, mConfig.OnlyDriverGenes ? mGeneDataCache : null)
                : mSampleVariantsCache.get(sampleId);

        if(refVariants.isEmpty())
            return;

        for(RefVariantData refVariant :refVariants)
        {
            // generate variant impact data and then write comparison results to CSV file
            VariantData variant = new VariantData(
                    refVariant.Chromosome, refVariant.Position, refVariant.Ref, refVariant.Alt);

            variant.setVariantDetails(refVariant.LocalPhaseSet, refVariant.Microhomology, refVariant.RepeatSequence, refVariant.RepeatCount);
            variant.setSampleId(sampleId);
            variant.setRefData(refVariant);

            try
            {
                variant.setRealignedVariant(createRightAlignedVariant(variant, mImpactClassifier.refGenome()));

                findVariantImpacts(variant, mImpactClassifier, mGeneDataCache);

                processPhasedVariants(variant.localPhaseSet());

                if(!variant.hasLocalPhaseSet())
                    processVariant(sampleId, variant, refVariant);
            }
            catch(Exception e)
            {
                PV_LOGGER.error("error processing var({})", variant);
                e.printStackTrace();
            }
        }

        processPhasedVariants(NO_LOCAL_PHASE_SET);
        mImpactClassifier.phasedVariants().clear();

        PV_LOGGER.debug("{}: sample({}) processed {} variants", mTaskId, sampleId, refVariants.size());
    }

    private void processPhasedVariants(int currentLocalPhaseSet)
    {
        List<VariantData> variants = mImpactClassifier.processPhasedVariants(currentLocalPhaseSet);

        if(variants != null)
            variants.forEach(x -> processVariant(x.sampleId(), x, x.refData()));
    }

    private void processVariant(final String sampleId, final VariantData variant, final RefVariantData refVariant)
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

        boolean hasDiff = false;

        if(variantImpact == null)
        {
            // hasDiff = refVariant.Reported;
            return;
        }
        else
        {
            if(mConfig.checkDiffType(REPORTED) && refVariant.Reported != variant.reported())
                hasDiff = true;

            if(mConfig.checkDiffType(CODING_EFFECT) && hasCodingEffectDiff(variantImpact, refVariant))
                hasDiff = true;

            VariantTransImpact transImpact = variant.getCanonicalTransImpacts(refVariant.Gene);

            if(transImpact != null)
            {
                VariantTransImpact raTransImpact = variant.getRealignedImpact(refVariant.Gene, transImpact);

                if(mConfig.checkDiffType(HGVS_CODING))
                {
                    if(hasHgvsCodingDiff(transImpact, raTransImpact, refVariant))
                        hasDiff = true;
                }

                if(mConfig.checkDiffType(HGVS_PROTEIN))
                {
                    if(hasHgvsProteinDiff(variant, transImpact, raTransImpact, refVariant))
                        hasDiff = true;
                }
            }
        }

        if(!hasDiff)
        {
            ++mMatchedCount;
        }
        else
        {
            logComparison(sampleId, refVariant, variant, variantImpact);
        }

        if(hasDiff || mConfig.WriteMatches)
            mWriter.writeVariantData(sampleId, variant, variantImpact, refVariant, hasDiff);
    }

    private boolean isReported(final VariantData variant, final VariantImpact variantImpact, final RefVariantData refVariant)
    {
        if(variantImpact == null)
            return false;

        if(mReportableOncoGenes.test(
                variantImpact.CanonicalGeneName, variant.type(), variant.repeatCount(), refVariant.IsHotspot,
                variantImpact.CanonicalCodingEffect, variantImpact.CanonicalEffect))
        {
            return true;
        }

        if(mReportableTsgGenes.test(
                variantImpact.CanonicalGeneName, variant.type(), variant.repeatCount(), refVariant.IsHotspot,
                variantImpact.CanonicalCodingEffect, variantImpact.CanonicalEffect))
        {
            return true;
        }

        return false;
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
                    sampleId, variant.toString(), variantImpact.CanonicalGeneName, variantImpact.CanonicalCodingEffect,
                    refVariant.Gene, refVariant.CanonicalCodingEffect);
            return;
        }

        PV_LOGGER.trace("sample({}) var({}) diff canonical coding: pave({} : {}) snpEff({} : {})",
                sampleId, variant.toString(), variantImpact.CanonicalGeneName, variantImpact.CanonicalCodingEffect,
                refVariant.Gene, refVariant.CanonicalCodingEffect);
    }

}
