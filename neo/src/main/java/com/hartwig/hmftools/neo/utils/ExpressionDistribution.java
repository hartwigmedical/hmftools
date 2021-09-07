package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.neo.utils.PeptideExpressionData.SOURCE_VALIDATION;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

public class ExpressionDistribution
{
    private final List<TpmBucket> mTpmBuckets;
    private final List<TpmBucket> mDecileBuckets;

    private final List<PeptideExpressionData> mAllExpData;

    private static final double MIN_TPM_BUCKET = 0.001;
    private static final double MAX_TPM_BUCKET = 100000;
    private static final double STRONG_BINDER_LIKELIHOOD = 0.001;

    public ExpressionDistribution(int totalPeptides)
    {
        mTpmBuckets = Lists.newArrayListWithCapacity(30);
        mTpmBuckets.add(new TpmBucket(0));

        double tpmBucket = MIN_TPM_BUCKET;
        while(tpmBucket < MAX_TPM_BUCKET)
        {
            mTpmBuckets.add(new TpmBucket(tpmBucket));
            tpmBucket *= 2;
        }

        mDecileBuckets = Lists.newArrayListWithCapacity(11);
        mAllExpData = Lists.newArrayListWithCapacity(totalPeptides);
    }

    public void process(final PeptideExpressionData pepExpData)
    {
        if(!pepExpData.hasTpm())
            return;

        if(pepExpData.Source.equals(SOURCE_VALIDATION) && pepExpData.LikelihoodRank > STRONG_BINDER_LIKELIHOOD)
            return;

        // add to ordered list for deciles and to TPM buckets
        addToTpmBucket(pepExpData);
        addToRankedList(pepExpData);
    }

    public void formTpmDeciles()
    {
        int peptidesPerDecile = mAllExpData.size() / 11;

        int nextCountTotal = peptidesPerDecile;

        TpmBucket currentBucket = new TpmBucket(0);
        mDecileBuckets.add(currentBucket);

        for(int i = 0; i < mAllExpData.size(); ++i)
        {
            PeptideExpressionData pepExpData = mAllExpData.get(i);

            if(i >= nextCountTotal || i == mAllExpData.size() - 1)
            {
                currentBucket.Bucket = pepExpData.tpm();

                if(mDecileBuckets.size() < 11)
                {
                    currentBucket = new TpmBucket(0);
                    mDecileBuckets.add(currentBucket);
                    nextCountTotal += peptidesPerDecile;
                }
            }

            currentBucket.add(pepExpData);
        }
    }

    private void addToRankedList(final PeptideExpressionData pepExpData)
    {
        // early exits
        if(mAllExpData.isEmpty())
        {
            mAllExpData.add(pepExpData);
            return;
        }

        double tpm = pepExpData.tpm();

        if(tpm < mAllExpData.get(0).tpm())
        {
            mAllExpData.add(0, pepExpData);
            return;
        }

        int itemCount = mAllExpData.size();
        if(tpm > mAllExpData.get(itemCount - 1).tpm())
        {
            mAllExpData.add(pepExpData);
            return;
        }

        if(itemCount < 20)
        {
            int index = 0;
            while(index < mAllExpData.size())
            {
                if(tpm < mAllExpData.get(index).tpm())
                    break;

                ++index;
            }

            mAllExpData.add(index, pepExpData);
            return;
        }

        int lowIndex = 0;
        int highIndex = mAllExpData.size() - 1;
        int currentIndex = mAllExpData.size() / 2;

        while(true)
        {
            double currentValue = mAllExpData.get(currentIndex).tpm();

            if(currentValue == tpm)
            {
                mAllExpData.add(currentIndex, pepExpData);
                return;
            }

            if(tpm < currentValue)
            {
                // current index is looking too high in the list
                if(currentIndex == lowIndex + 1)
                {
                    // no need to look any lower (again
                    mAllExpData.add(currentIndex, pepExpData);
                    return;
                }

                highIndex = currentIndex;
            }
            else
            {
                if(currentIndex == highIndex - 1)
                {
                    mAllExpData.add(currentIndex + 1, pepExpData);
                    return;
                }

                lowIndex = currentIndex;
            }

            int newIndex = lowIndex + (highIndex - lowIndex) / 2;

            if(newIndex == currentIndex)
            {
                mAllExpData.add(currentIndex, pepExpData);
                return;
            }

            currentIndex = newIndex;
        }
    }

    private void addToTpmBucket(final PeptideExpressionData pepExpData)
    {
        double tpm = pepExpData.tpm();

        for(int i = 0; i < mTpmBuckets.size(); ++i)
        {
            TpmBucket bucket = mTpmBuckets.get(i);

            if(tpm > bucket.Bucket)
                continue;

            if(tpm < bucket.Bucket || Doubles.equal(bucket.Bucket, tpm))
            {
                bucket.add(pepExpData);
                break;
            }

            if(i < mTpmBuckets.size() - 1)
            {
                TpmBucket nextBucket = mTpmBuckets.get(i + 1);

                if(tpm < nextBucket.Bucket || Doubles.equal(nextBucket.Bucket, tpm))
                {
                    nextBucket.add(pepExpData);
                    break;
                }
            }
            else
            {
                TpmBucket lastBucket = mTpmBuckets.get(mTpmBuckets.size() - 1);
                lastBucket.add(pepExpData);
            }
        }
    }

    private class TpmBucket
    {
        public double Bucket;
        public int TotalCount;
        public double StrongBinders;

        public TpmBucket(double bucket)
        {
            Bucket = bucket;
            TotalCount = 0;
            StrongBinders = 0;
        }

        public void add(final PeptideExpressionData pepExpData)
        {
            ++TotalCount;

            if(pepExpData.Source.equals(SOURCE_VALIDATION))
                ++StrongBinders;
        }

        public double calcRate() { return TotalCount > 0 ? StrongBinders / (double)TotalCount : 0; }
    }
}
