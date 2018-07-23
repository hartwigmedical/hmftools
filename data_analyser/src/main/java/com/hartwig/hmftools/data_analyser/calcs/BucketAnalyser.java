package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.round;

import static com.hartwig.hmftools.data_analyser.calcs.NmfConfig.NMF_FS_MIN_SAMPLES;
import static com.hartwig.hmftools.data_analyser.types.BucketPair.RATIO_MAX;
import static com.hartwig.hmftools.data_analyser.types.BucketPair.RATIO_MIN;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.types.BucketGroup;
import com.hartwig.hmftools.data_analyser.types.BucketPair;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BucketAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(NmfManager.class);

    private String mOutputDir;

    private GenericDataCollection mDataCollection;

    private NmfMatrix mSampleCountsMatrix;

    // a matrix with results from comparing pairs of buckets, with pairingId down the rows and sample index in the cols
    private NmfMatrix mBPMatrix;

    // a link between the BP id and the 2 bucket indices it relates to
    private List<int[]> mBPIds;

    private List<BucketPair> mBucketPairs;
    private List<BucketGroup> mBucketGroups;
    private int mMinSampleCount;

    private static int MIN_BG_BUCKET_COUNT = 20; // below which a bucket ratio is considered invalid

    public BucketAnalyser()
    {
        mOutputDir = "";
        mDataCollection = null;
        mSampleCountsMatrix = null;
        mMinSampleCount = 0;

        mBPMatrix = null;
        mBPIds = null;
        mBucketPairs = Lists.newArrayList();
        mBucketGroups = Lists.newArrayList();
    }

    public void setOutputDir(final String outputDir)
    {
        mOutputDir = outputDir;
    }

    public static void addCmdLineArgs(Options options) {

    }

    public void initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        // mConfig = new NmfConfig(cmd);
        mDataCollection = collection;
        // bucketAnalyser.setOutputDir(cmd.getOptionValue(OUTPUT_DIR));

        double minSamplePerc = cmd.hasOption(NMF_FS_MIN_SAMPLES) ? Double.parseDouble(cmd.getOptionValue(NMF_FS_MIN_SAMPLES)) : 0.01;

        mSampleCountsMatrix = DataUtils.createMatrixFromListData(mDataCollection.getData());
        mSampleCountsMatrix.cacheTranspose();

        mMinSampleCount = (int)round(minSamplePerc * mSampleCountsMatrix.Cols);

        LOGGER.info("bucketCount({}) sampleCount({})", mSampleCountsMatrix.Rows, mSampleCountsMatrix.Cols);

        // report bucket info
        for(int i = 0; i < mSampleCountsMatrix.Rows; ++i)
        {
            int bucketTotal = (int)DataUtils.sumVector(mSampleCountsMatrix.getRow(i));
            LOGGER.debug("bucket({}) count({})", i, bucketTotal);
        }

    }

    public void run()
    {
        PerformanceCounter perfCounter = new PerformanceCounter("BucketAnalysis");

        perfCounter.start("Ratios");
        calcBucketRatios();
        perfCounter.stop();

        perfCounter.start("Pairings");
        createBucketPairs();
        perfCounter.stop();

        perfCounter.start("LinkPairs");
        linkBucketPairs();
        perfCounter.stop();

        perfCounter.start("Groups");
        formBucketGroups();
        perfCounter.stop();

        perfCounter.logStats(true);
    }

    private void calcBucketRatios()
    {
        int bucketCount = mSampleCountsMatrix.Rows;
        int sampleCount = mSampleCountsMatrix.Cols;
        int pairingIndex = 0;

        int bucketPairingCount = bucketCount  * (bucketCount - 1) / 2;
        mBPMatrix = new NmfMatrix(bucketPairingCount, sampleCount);
        mBPIds = Lists.newArrayList();
        double[][] pairingsData = mBPMatrix.getData();
        final double[][] bucketCountData = mSampleCountsMatrix.getData();

        LOGGER.info("calculating ratios for {} buckets, {} samples, total = {}",
                mSampleCountsMatrix.Rows, mSampleCountsMatrix.Cols, bucketPairingCount);

        for (int a = 0; a < bucketCount; ++a) {

            for (int b = a + 1; b < bucketCount; ++b) {

                for(int n = 0; n < sampleCount; ++n) {

                    double countA = bucketCountData[a][n];
                    double countB = bucketCountData[b][n];
                    double bpCountRatio = calcRatio(countA, countB);

                    pairingsData[pairingIndex][n] = bpCountRatio;
                }

                int[] pairingId = {a,b};
                mBPIds.add(pairingId);
                ++pairingIndex;
            }
        }
    }

    private double calcRatio(double countA, double countB)
    {
        if(countA <= 0 || countB <= 0)
            return 0;

        // consider a ratio unreliable for smaller bucket counts
        if(countA <= MIN_BG_BUCKET_COUNT || countB <= MIN_BG_BUCKET_COUNT)
            return 0;

        // ratio is lower index bucket over higher
        return Math.max(Math.min(countA / countB, RATIO_MAX), RATIO_MIN);
    }

    private void createBucketPairs()
    {
        // look for repeated or similar ratios across the samples for the same bucket pairings
        // keep track of frequencies of BP ratios

        // List<HashMap<Double,Integer>> mBPFrequencies

        final double[][] pairingsData = mBPMatrix.getData();
        int sampleBucketArea = 0;

        for(int bpIndex = 0; bpIndex < mBPMatrix.Rows; ++bpIndex)
        {
            // for this bucket pairing, check all samples for similarities in their ratios,
            // and where found, add those sampleIds to the same BucketPairing
            List<BucketPair> bucketPairs = Lists.newArrayList();

            for(int sampleIndex = 0; sampleIndex < mBPMatrix.Cols; ++sampleIndex)
            {
                double bpRatio = pairingsData[bpIndex][sampleIndex];

                if(bpRatio <= RATIO_MIN || bpRatio >= RATIO_MAX)
                    continue;

                boolean found = false;
                for(BucketPair bucketPair : bucketPairs)
                {
                    if(bucketPair.isWithinRange(bpRatio))
                    {
                        bucketPair.addSampleRatio(sampleIndex, bpRatio);
                        found = true;
                        break;
                    }
                }

                if(!found)
                {
                    int[] bucketIds = mBPIds.get(bpIndex);

                    BucketPair bucketPair = new BucketPair(bpIndex, bucketIds[0], bucketIds[1]);
                    bucketPair.addSampleRatio(sampleIndex, bpRatio);
                    bucketPairs.add(bucketPair);
                }
            }

            int samplesPaired = 0;
            for(BucketPair bucketPair : bucketPairs)
            {
                if(bucketPair.getSampleIds().size() >= mMinSampleCount) {

                    mBucketPairs.add(bucketPair);
                    sampleBucketArea += bucketPair.getSampleIds().size();
                    samplesPaired += bucketPair.getSampleIds().size();

                    LOGGER.debug(String.format("added pair(%d [%d][%d]: samples(%d) ratio(%.3f range=%.3f -> %.3f asPerc=%.3f) bucketPairCount(%d) areaPaired(%d)",
                            bucketPair.getId(), bucketPair.getBucketA(), bucketPair.getBucketB(), bucketPair.getSampleIds().size(),
                            bucketPair.getAverageRatio(), bucketPair.getLowRatio(), bucketPair.getHighRatio(), bucketPair.calcRangePercent(),
                            mBucketPairs.size(), sampleBucketArea));
                }
            }

            if(samplesPaired > 0) {
                double samplePairedPerc = samplesPaired / (double) mBPMatrix.Cols;
                LOGGER.debug(String.format("bucketPair(%d) totalSamples(%d perc=%.3f)", bpIndex, samplesPaired, samplePairedPerc));
            }
        }

        int totalArea = mBPMatrix.Rows * mBPMatrix.Cols;
        double totalPairedPerc = sampleBucketArea / (double)totalArea;

        LOGGER.debug(String.format("total areaPaired(%d) of total(%d) percent(%.3f)",
                sampleBucketArea, totalArea, totalPairedPerc));
    }

    private void linkBucketPairs()
    {
        // merge any bucket pairs which have over-lapping ratios
        int initBucketPairs = mBucketPairs.size();

        int index1 = 0;

        while(index1 < mBucketPairs.size())
        {
            BucketPair bp1 = mBucketPairs.get(index1);

            // check in the remaining set of bucket pairs
            for(int index2 = index1 + 1; index2 < mBucketPairs.size(); ++index2) {

                final BucketPair bp2 = mBucketPairs.get(index2);

                // must be the same pair
                if (bp1.getId() != bp2.getId())
                {
                    continue;
                }

                // ovelapping ratios
                if(bp1.getHighRatio() <= bp2.getLowRatio() || bp1.getLowRatio() >= bp2.getHighRatio())
                {
                    continue;
                }

                // join them up
                LOGGER.debug(String.format("joining pairs(%d: %d & %d) samples(%d & %d) ratios(%.3f -> %.3f & %.3f -> %.3f)",
                        bp1.getId(), bp1.getBucketA(), bp1.getBucketB(), bp1.getSampleIds().size(), bp2.getSampleIds().size(),
                        bp1.getLowRatio(), bp1.getHighRatio(), bp2.getLowRatio(), bp2.getHighRatio()));

                for(int j = 0; j < bp2.getSampleIds().size(); ++j) {

                    bp1.addSampleRatio(bp2.getSampleIds().get(j), bp2.getRatios().get(j));
                }

                LOGGER.debug(String.format("composite pair(%d [%d][%d]: samples(%d) ratio(%.3f range=%.3f -> %.3f asPerc=%.3f)",
                        bp1.getId(), bp1.getBucketA(), bp1.getBucketB(), bp1.getSampleIds().size(),
                        bp1.getAverageRatio(), bp1.getLowRatio(), bp1.getHighRatio(), bp1.calcRangePercent()));

                mBucketPairs.remove(index2);
                break;
            }

            ++index1;

        } // end for each BP using the first index

        LOGGER.info("bucket pairs count after similar ratio joining({} -> {})", initBucketPairs, mBucketPairs.size());
    }

    private void formBucketGroups()
    {
        if(mBucketPairs.isEmpty())
            return;

        LOGGER.debug("creating groups from {} bucket pairs", mBucketPairs.size());

        // look for linked bucket ids across groups containing the same subset of samples

        int index1 = 0;
        int index2 = 0;
        while(index1 < mBucketPairs.size())
        {
            if(index1 > 0 && (index1 % 100) == 0)
            {
                LOGGER.debug("processed {} of {} bucket pairs", index1, mBucketPairs.size());
            }

            final BucketPair bp1 = mBucketPairs.get(index1);

            boolean pairAdded = false;

            // first check existing groups
            for(BucketGroup bucketGroup : mBucketGroups)
            {
                // check for overlapping buckets
                if (!bucketGroup.hasBucket(bp1.getBucketA()) && !bucketGroup.hasBucket(bp1.getBucketB()))
                    continue;

                // and then overlapping samples
                List<Integer> sharedSamples = bp1.getSharedSamples(bucketGroup.getSampleIds());
                int sharedCount = sharedSamples.size();

                if(sharedCount < mMinSampleCount)
                    continue;

                // add this BP to existing group
                bucketGroup.addBucketPair(bp1);
                // bucketGroup.addSamples(sharedSamples); // no new samples, since only taken the union with existing

                LOGGER.debug("bucketGroup({}) added pair({} shared({} of {}) samples({}) buckets({}) pairs({})",
                        bucketGroup.getId(), bp1.getId(), sharedCount, bp1.getSampleIds().size(),
                        bucketGroup.getSampleIds().size(), bucketGroup.getBucketIds().size(),
                        bucketGroup.getBucketPairs().size());

                pairAdded = true;
                break;
            }

            if(pairAdded)
            {
                mBucketPairs.remove(index1);
                continue;
            }

            // and if not matched, check in the remaining set of bucket pairs
            index2 = index1 + 1;
            while(index2 < mBucketPairs.size()) {

                final BucketPair bp2 = mBucketPairs.get(index2);

                // check that these bucket pairs aren't for exactly the same bucket-pair
                // (which would implying they have diff ratios, and so different sets of samples)
                if (bp1.getId() == bp2.getId())
                {
                    ++index2;
                    continue;
                }

                // check for overlapping buckets
                if (!bp1.hasSameBucket(bp2))
                {
                    ++index2;
                    continue;
                }

                // and then overlapping samples
                List<Integer> sharedSamples = bp1.getSharedSamples(bp2.getSampleIds());

                int sharedCount = sharedSamples.size();

                if (sharedCount < mMinSampleCount)
                {
                    ++index2;
                    continue;
                }

                // form a new bucket group
                BucketGroup bucketGroup = new BucketGroup(mBucketGroups.size(), sharedSamples);
                bucketGroup.addBucketPair(bp1);
                bucketGroup.addBucketPair(bp2);
                mBucketGroups.add(bucketGroup);

                LOGGER.debug("added new bucketGroup({}) with pairs({} & {}) samples(shared={} of bp1={} & bp2={})",
                        bucketGroup.getId(), index1, index2, bucketGroup.getSampleIds().size(),
                        bp1.getSampleIds().size(), bp2.getSampleIds().size());

                pairAdded = true;
                break;
            }

            if(pairAdded)
            {
                mBucketPairs.remove(index2); // need to remove later index first
                mBucketPairs.remove(index1);
            }
            else
            {
                // check next item for the first index run
                ++index1;
            }

        } // end for each BP using the first index

        // report the final set of buckets
        LOGGER.info("summary of {} bucketGroups:", mBucketGroups.size());

        int areaGrouped = 0;

        for (BucketGroup bucketGroup : mBucketGroups) {

            int sampleCount = bucketGroup.getSampleIds().size();
            int bucketCount = bucketGroup.getBucketIds().size();

            LOGGER.info("bucketGroup({}) samples({}) buckets({}: {}) pairs({})",
                    bucketGroup.getId(), sampleCount, bucketCount, bucketGroup.getBucketIds().toString(),
                    bucketGroup.getBucketPairs().size());

            areaGrouped += sampleCount * bucketCount;
        }

        double totalArea = mSampleCountsMatrix.Cols * mSampleCountsMatrix.Rows;

        LOGGER.info(String.format("areaGrouped(%d) total(%.0f) percent(%.3f)",
                areaGrouped, totalArea, areaGrouped/totalArea));
    }
}
