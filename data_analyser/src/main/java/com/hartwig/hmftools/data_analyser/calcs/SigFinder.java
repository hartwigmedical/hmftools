package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I1;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I2;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_VAL;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.calcCSS;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getTopCssPairs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;
import com.hartwig.hmftools.data_analyser.types.SigGroup;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigFinder {

    final NmfMatrix mSampleCounts;
    final NmfConfig mConfig;
    double mTotalVariants;

    List<SigGroup> mSigGroups;
    NmfMatrix mSignatures;
    NmfMatrix mContributions;
    NmfMatrix mAllocatedCounts; // a working contribution matrix for assiging counts from proposed sigs
    List<Integer> mProposedSigs;

    final String mOutputDir;
    final String mOutputFileId;
    final List<String> mSampleNames;

    private static final Logger LOGGER = LogManager.getLogger(SigFinder.class);

    public SigFinder(
            final NmfMatrix sampleBucketCounts, final NmfConfig config,
            final String outputDir,final String outputFileId, final List<String> sampleNames)
    {
        mSampleCounts = sampleBucketCounts;
        mConfig = config;
        mOutputDir = outputDir;
        mOutputFileId = outputFileId;
        mSampleNames = sampleNames;
        mProposedSigs = Lists.newArrayList();

        mTotalVariants = mSampleCounts.sum();

        mSigGroups = Lists.newArrayList();
        mSignatures = null;
        mContributions = null;
        mAllocatedCounts = null;
    }

    public final NmfMatrix getSignatures() { return mSignatures; }
    public final NmfMatrix getContributions() { return mContributions; }

    public void findSignatures()
    {
        LOGGER.debug("finding signatures");

        PerformanceCounter perfCounter = new PerformanceCounter("SigFinder");

        perfCounter.start();
        formSigGroups();
        proposeSignatures();
        // applyProposedSignatures();
        createSignatures();

        perfCounter.stop();

        if(mSignatures == null)
            return;

        CosineSim.logSimilarites(mSignatures, 0.8, "sig");

        if(mSignatures != null && mContributions != null
        && mOutputDir != null && !mSampleNames.isEmpty())
        {
            writeSignatures(mSignatures);
            writeContributions(mContributions);
        }

        perfCounter.logStats();
    }

    private double calcLogLikelihood(int sam1, int sam2)
    {
        double[] sData1 = mSampleCounts.getCol(sam1);
        double[] sData2 = mSampleCounts.getCol(sam2);

        return CosineSim.calcLogLikelihood(sData1, sData2, true);
    }

    private void formSigGroups()
    {
        // run CSS on all pairs, storing the highest value ones (not mutually exclusive yet)
        int bucketCount = mSampleCounts.Rows;

        final List<double[]> cssResults = getTopCssPairs(mSampleCounts, mSampleCounts, mConfig.CssCutoff, false, true);

        // final List<double[]> lliResults = getTopLogLikelihoodPairs(mSampleCounts, mSampleCounts, mConfig.CssCutoff);

        // create potential signature groups, only allocating a sample to at most 1 group
        List<Integer> allocatedSampleIds = Lists.newArrayList();

        for(final double[] cssResult : cssResults)
        {
            int samId1 = (int)cssResult[CSSR_I1];
            int samId2 = (int)cssResult[CSSR_I2];

            if(allocatedSampleIds.contains(samId1) && allocatedSampleIds.contains(samId2))
                continue;

            // compare with log-likelihood probability
            // double lliProb = calcLogLikelihood(samId1, samId2);

            // check for a group with either of these samples in it already
            boolean sampleFound = false;
            for(SigGroup sigGroup : mSigGroups)
            {
                boolean hasSample1 = sigGroup.hasSampleId(samId1);
                boolean hasSample2 = sigGroup.hasSampleId(samId2);

                if(!hasSample1 && !hasSample2 || (hasSample1 && hasSample2))
                    continue;

                sampleFound = true;

                // calculate and test a new CSS vs the existing group
                int newSampleId = !hasSample1 ? samId1 : samId2;

                double[] sampleCounts = mSampleCounts.getCol(newSampleId);
                double newCss = calcCSS(sampleCounts, sigGroup.getBucketCounts());

                if(newCss >= mConfig.CssCutoff)
                {
                    sigGroup.addSample(newSampleId, sampleCounts, min(newCss, cssResult[CSSR_VAL]));
                    allocatedSampleIds.add(newSampleId);

                    LOGGER.debug(String.format("sigGroup(%d) adding sample(%d) groupSamples(%d) css(%.4f vsGrp=%.4f init=%.4f) totalCount(%d)",
                            sigGroup.id(), newSampleId, sigGroup.getSampleIds().size(),
                            cssResult[CSSR_VAL], newCss, sigGroup.getInitialCss(), sigGroup.getTotalCount()));
                }
                else
                {
                    LOGGER.debug(String.format("sigGroup(%d) skipping sample(%d) css(%.4f vsGrp=%.4f init=%.4f)",
                            sigGroup.id(), newSampleId, cssResult[CSSR_VAL], newCss, sigGroup.getInitialCss()));
                }

                break;
            }

            // if one or the other sample is already part of a group but the unallocated sample wasn't a close enough match,
            // continue and don't create a new group with on the already allocated samples
            if(sampleFound)
                continue;

            SigGroup sigGroup = new SigGroup(mSigGroups.size(), bucketCount, cssResult[CSSR_VAL]);
            sigGroup.addSample(samId1, mSampleCounts.getCol(samId1), cssResult[CSSR_VAL]);
            sigGroup.addSample(samId2, mSampleCounts.getCol(samId2), cssResult[CSSR_VAL]);
            allocatedSampleIds.add(samId1);
            allocatedSampleIds.add(samId2);

            LOGGER.debug(String.format("adding sigGroup(%d) samples(%d & %d) css(%.4f) totalCount(%d)",
                    sigGroup.id(), samId1, samId2, cssResult[CSSR_VAL], sigGroup.getTotalCount()));

            mSigGroups.add(sigGroup);
        }

        // reconcile groups
        int index1 = 0;
        while(index1 < mSigGroups.size())
        {
            SigGroup sg1 = mSigGroups.get(index1);

            int index2 = index1 + 1;
            while(index2 < mSigGroups.size())
            {
                SigGroup sg2 = mSigGroups.get(index2);

                double newCss = calcCSS(sg1.getBucketCounts(), sg2.getBucketCounts());

                if(newCss >= mConfig.CssCutoff)
                {
                    double minCss = min(newCss, sg2.getWorstCss());

                    for(Integer sampleId : sg2.getSampleIds()) {
                        sg1.addSample(sampleId, mSampleCounts.getCol(sampleId), minCss);
                    }

                    LOGGER.debug(String.format("sigGroups(%d & %d) joined groupSamples(%d) css(%.4f init=%.4f) totalCount(%d)",
                            sg1.id(), sg2.id(), sg1.getSampleIds().size(), newCss, sg1.getInitialCss(), sg1.getTotalCount()));

                    mSigGroups.remove(index2);
                    continue;
                }

                ++index2;
            }

            ++index1;
        }

        if(mSigGroups.isEmpty())
        {
            LOGGER.info("no sig groups found");
            return;
        }

        LOGGER.info("created {} sigGroups allocating {} samples",
                mSigGroups.size(), allocatedSampleIds.size());

        for(final SigGroup sigGroup : mSigGroups)
        {
            double groupPercent = sigGroup.getTotalCount() / mTotalVariants;

            LOGGER.debug(String.format("sigGroup(%d) sampleCount(%d) totalCount(%d perc=%.2f) css(init=%.4f worst=%.4f)",
                    sigGroup.id(), sigGroup.getSampleIds().size(), sigGroup.getTotalCount(), groupPercent,
                    sigGroup.getInitialCss(), sigGroup.getWorstCss()));
        }
    }

    private void proposeSignatures()
    {
        if (mSigGroups.isEmpty())
            return;

        int totalCountIncluded = 0;
        int totalSamplesIncluded = 0;

        int allSampleCount = mSampleCounts.Cols;
        double minSampleCount = mConfig.MinSamplePerc * allSampleCount;
        List<Integer> tmpProposedSigs = Lists.newArrayList();

        List<Double> proposedSigScores = Lists.newArrayList();

        for (int i = 0; i < mSigGroups.size(); ++i)
        {
            final SigGroup sigGroup = mSigGroups.get(i);

            // skip similarities only involve a few samples
            int groupSampleCount = sigGroup.getSampleIds().size();

            if (groupSampleCount < minSampleCount)
                continue;

            // skip below significant similarity in samples
            double worstCss = sigGroup.getWorstCss();

            if (worstCss < mConfig.CssCutoff)
                continue;

            double groupVarPerc = sigGroup.getTotalCount() / mTotalVariants;

            // skip if this signature will only affect a small percent of variants
            if (groupVarPerc < mConfig.MinSamplePerc)
                continue;

            double groupSamplePerc = groupSampleCount / (double)allSampleCount;
            double cssFactor = (worstCss - mConfig.CssCutoff) / (1 - mConfig.CssCutoff);
            double sigScore = cssFactor * groupSamplePerc * groupVarPerc;
            sigGroup.setScore(sigScore);

            tmpProposedSigs.add(i);
            proposedSigScores.add(sigScore);
            totalCountIncluded += sigGroup.getTotalCount();
            totalSamplesIncluded += sigGroup.getSampleIds().size();

            LOGGER.info(String.format("sigGroup(%d) proposed with score(%.4g) sampleCount(%d perc=%.3f) totalCount(%d perc=%.3f) css(init=%.4f worst=%.4f)",
                    sigGroup.id(), sigScore, groupSampleCount, groupSamplePerc, sigGroup.getTotalCount(),
                    groupVarPerc, sigGroup.getInitialCss(), worstCss));

        }

        int sigCount = tmpProposedSigs.size();

        if (sigCount == 0)
        {
            LOGGER.info("no valid proposed signatures found");
            return;
        }
        else if (sigCount > mConfig.SigCount)
        {
            LOGGER.warn("more proposed sigs({}) than permitted({})", sigCount, mConfig.SigCount);

            // order the proposed sigs by score and take the top X
            double[] psValues = new double[proposedSigScores.size()];
            for(int i = 0; i < psValues.length; ++i)
                psValues[i] = proposedSigScores.get(i);

            List<Integer> sortedScoredIndices = getSortedVectorIndices(psValues, false);

            for(int i = 0; i < mConfig.SigCount; ++i)
            {
                int proposedSigIndex = sortedScoredIndices.get(i);
                mProposedSigs.add(tmpProposedSigs.get(proposedSigIndex));
            }
        }
        else
        {
            mProposedSigs = tmpProposedSigs;
        }

        double totalCountPercent = totalCountIncluded / mTotalVariants;
        double totalSamplesPercent = totalSamplesIncluded / (double) mSampleCounts.Cols;

        LOGGER.info(String.format("%d proposed sigs, totalCount(%d asPerc=%.3f) totalSamples(%d asPerc=%.3f)",
                sigCount, totalCountIncluded, totalCountPercent, totalSamplesIncluded, totalSamplesPercent));
    }

    private void applyProposedSignatures()
    {
        int bucketCount = mSampleCounts.Rows;
        int sampleCount = mSampleCounts.Cols;

        mAllocatedCounts = new NmfMatrix(bucketCount, sampleCount);
        double[][] aData = mAllocatedCounts.getData();

        for(Integer sigGroupId : mProposedSigs)
        {
            final SigGroup sigGroup = mSigGroups.get(sigGroupId);
            final double[] bucketRatios = sigGroup.getBucketRatios();

            for(int n = 0; n < sampleCount; ++n)
            {
                if(!sigGroup.hasSampleId(n))
                    continue; // value left as 1

                double lowestBucketFactor = 0;

                double[] sampleCounts = mSampleCounts.getCol(n);

                for(int i = 0; i < bucketCount; ++i)
                {
                    if(bucketRatios[i] == 0)
                        continue;

                    // a bucket count of zero is bumped up to 1 on the assumption that
                    double sbRatio = max(sampleCounts[i],1) / bucketRatios[i];

                    if(sbRatio > 0 && (i == 0 || sbRatio < lowestBucketFactor))
                        lowestBucketFactor = sbRatio;
                }

                // now allocate/remove the corresponding counts from this sample
                for(int i = 0; i < bucketCount; ++i)
                {
                    double allocatedCount = floor(lowestBucketFactor * bucketRatios[i]);

                    if(allocatedCount == 0)
                        continue;

                    aData[i][n] += allocatedCount;

                    LOGGER.debug(String.format("sample(%d) bucket(%d) allocated(%.0f of %.0f)",
                            n, i, allocatedCount, sampleCounts[i]));
                }
            }
        }
    }

    private void createSignatures()
    {
        if(mProposedSigs.isEmpty())
            return;

        int bucketCount = mSampleCounts.Rows;
        int allSampleCount = mSampleCounts.Cols;
        int sigCount = mProposedSigs.size();

        mSignatures = new NmfMatrix(bucketCount, sigCount);
        double[][] sigData = mSignatures.getData();

        mContributions = new NmfMatrix(sigCount, mSampleCounts.Cols);
        double[][] contribData = mContributions.getData();

        int sigIndex = 0;
        for(final Integer sigGroupId : mProposedSigs)
        {
            final SigGroup sigGroup = mSigGroups.get(sigGroupId);

            LOGGER.info(String.format("proposed sig(%d) score(%.4g) sampleCount(%d perc=%.3f) totalCount(%d perc=%.3f) css(%.4f)",
                    sigIndex, sigGroup.getScore(), sigGroup.getSampleIds().size(), sigGroup.getSampleIds().size() / (double)allSampleCount,
                    sigGroup.getTotalCount(), sigGroup.getTotalCount() / mTotalVariants, sigGroup.getWorstCss()));

            final double[] bucketRatios = sigGroup.getBucketRatios();

            for(int j = 0; j < bucketCount; ++j)
            {
                sigData[j][sigIndex] = bucketRatios[j];
            }

            for(int n = 0; n < mContributions.Cols; ++n)
            {
                if(!sigGroup.hasSampleId(n))
                    continue; // value left as 1

                // check CSS again in case it's deviated from this sample's counts too much
                double css = calcCSS(sigGroup.getBucketCounts(), mSampleCounts.getCol(n));

                if(css < mConfig.CssCutoff)
                    continue;

                contribData[sigIndex][n] = 1;
            }

            // log the top 10 contributions of any significance
            final double[] sigValues = mSignatures.getCol(sigIndex);
            List<Integer> sortedIndices = DataUtils.getSortedVectorIndices(sigValues, false);

            for(int i = 0; i < 5; ++i)
            {
                int bucketIndex = sortedIndices.get(i);
                LOGGER.debug(String.format("proposed sig(%d) bucket: %d = %.3f", sigIndex, bucketIndex, sigValues[bucketIndex]));

                if(sigValues[bucketIndex] < 0.01)
                    break;
            }

            ++sigIndex;
        }
    }

    public void writeSignatures(final NmfMatrix signatures)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir,mOutputFileId + "_sf_sigs.csv");

            int i = 0;
            for(; i < signatures.Cols-1; ++i)
            {
                writer.write(String.format("%d,", i));
            }
            writer.write(String.format("%d", i));

            writer.newLine();

            writeMatrixData(writer, signatures, false);

            writer.close();
        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void writeContributions(final NmfMatrix contributions)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_sf_contribs.csv");

            int i = 0;
            for(; i < mSampleNames.size()-1; ++i)
            {
                writer.write(String.format("%s,", mSampleNames.get(i)));
            }
            writer.write(String.format("%s", mSampleNames.get(i)));

            writer.newLine();

            writeMatrixData(writer, contributions, false);

            writer.close();
        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }


}
