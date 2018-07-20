package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.DataAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I1;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I2;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_VAL;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.calcCSS;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getTopCssPairs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import javax.xml.crypto.Data;

import com.google.common.collect.Lists;
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

        mTotalVariants = mSampleCounts.sum();

        mSigGroups = Lists.newArrayList();
        mSignatures = null;
        mContributions = null;
    }

    public final NmfMatrix getSignatures() { return mSignatures; }
    public final NmfMatrix getContributions() { return mContributions; }

    public void findSignatures()
    {
        LOGGER.debug("finding signatures");

        formSigGroups();
        createSignatures();

        if(mSignatures != null)
            CosineSim.logSimilarites(mSignatures, 0.8, "sig");

        if(mSignatures != null && mContributions != null
        && mOutputDir != null && !mSampleNames.isEmpty())
        {
            writeSignatures(mSignatures);
            writeContributions(mContributions);
        }
    }

    private void formSigGroups()
    {
        // run CSS on all pairs, storing the highest value ones (not mutually exclusive yet)
        int bucketCount = mSampleCounts.Rows;

        final List<double[]> cssResults = getTopCssPairs(mSampleCounts, mSampleCounts, mConfig.CssCutoff, false, true);

        // create potential signature groups, only allocating a sample to at most 1 group
        List<Integer> allocatedSampleIds = Lists.newArrayList();

        for(final double[] cssResult : cssResults)
        {
            int samId1 = (int)cssResult[CSSR_I1];
            int samId2 = (int)cssResult[CSSR_I2];

            if(allocatedSampleIds.contains(samId1) && allocatedSampleIds.contains(samId2))
                continue;

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

        LOGGER.info("created {} sigGroups allocating {} samples",
                mSigGroups.size(), allocatedSampleIds.size());

        for(final SigGroup sigGroup : mSigGroups)
        {
            double groupPercent = sigGroup.getTotalCount() / mTotalVariants;

            LOGGER.debug(String.format("sigGroup(%d) sampleCount(%d) totalCount(%d asPerc=%.2f) css(init=%.4f worst=%.4f)",
                    sigGroup.id(), sigGroup.getSampleIds().size(), sigGroup.getTotalCount(), groupPercent,
                    sigGroup.getInitialCss(), sigGroup.getWorstCss()));
        }
    }

    private void createSignatures()
    {
        if(mSigGroups.isEmpty())
            return;

        int bucketCount = mSampleCounts.Rows;
        int totalCountIncluded = 0;
        int totalSamplesIncluded = 0;

        List<Integer> proposedSigs = Lists.newArrayList();

        double percentCutoff = 0.03; // of the total variants

        for(int i = 0; i < mSigGroups.size(); ++i)
        {
            final SigGroup sigGroup = mSigGroups.get(i);

            // skip similarities only involve a few samples
            if(sigGroup.getSampleIds().size() < mConfig.MinSampleCount)
                continue;

            // skip below significant similarity in samples
            if(sigGroup.getWorstCss() < mConfig.CssCutoff)
                continue;

            double groupPercent = sigGroup.getTotalCount() / mTotalVariants;

            // skip if this signature will only affect a small percent of variants
            if(groupPercent < percentCutoff)
                continue;

            proposedSigs.add(i);
            totalCountIncluded += sigGroup.getTotalCount();
            totalSamplesIncluded += sigGroup.getSampleIds().size();

            LOGGER.info(String.format("sigGroup(%d) proposed with sampleCount(%d) totalCount(%d asPerc=%.3f) css(init=%.4f worst=%.4f)",
                    sigGroup.id(), sigGroup.getSampleIds().size(), sigGroup.getTotalCount(),
                    groupPercent, sigGroup.getInitialCss(), sigGroup.getWorstCss()));

            if(proposedSigs.size() >= mConfig.SigCount)
                break;
        }

        int sigCount = proposedSigs.size();

        if(sigCount == 0)
        {
            LOGGER.info("no valid proposed signatures found");
            return;
        }
        else if(sigCount > mConfig.SigCount)
        {
            LOGGER.warn("more proposed sigs({}) than permitted({}), need a ranking system", sigCount, mConfig.SigCount);
        }

        double totalCountPercent = totalCountIncluded / mTotalVariants;
        double totalSamplesPercent = totalSamplesIncluded / (double)mSampleCounts.Cols;

        LOGGER.info(String.format("%d proposed sigs, totalCount(%d asPerc=%.3f) totalSamples(%d asPerc=%.3f)",
                sigCount, totalCountIncluded, totalCountPercent, totalSamplesIncluded, totalSamplesPercent));

        mSignatures = new NmfMatrix(bucketCount, sigCount);
        double[][] sigData = mSignatures.getData();

        mContributions = new NmfMatrix(sigCount, mSampleCounts.Cols);
        double[][] contribData = mContributions.getData();

        int sigIndex = 0;
        for(final Integer sigGroupId : proposedSigs)
        {
            final SigGroup sigGroup = mSigGroups.get(sigGroupId);
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
            List<Integer> sortedList = DataUtils.getSortedVector(sigValues, false);

            for(int i = 0; i < 5; ++i)
            {
                int bucketIndex = sortedList.get(i);
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
