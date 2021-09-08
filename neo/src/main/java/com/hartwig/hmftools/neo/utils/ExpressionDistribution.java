package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.EXP_TYPE_DECILE;
import static com.hartwig.hmftools.neo.bind.BindCommon.EXP_TYPE_TPM_LEVEL;
import static com.hartwig.hmftools.neo.bind.BindCommon.formFilename;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_EXPRESSION_DIST;
import static com.hartwig.hmftools.neo.utils.PeptideExpressionData.SOURCE_VALIDATION;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

public class ExpressionDistribution
{
    private final List<TpmBucket> mTpmBuckets;
    private final List<TpmBucket> mDecileBuckets;

    private final List<TpmBucket> mRoundedTpmBuckets;

    private static final int MIN_TPM_EXP = -5;
    private static final int MAX_TPM_EXP = 15;
    private static final double MIN_TPM_BUCKET = pow(2, MIN_TPM_EXP);
    private static final double MAX_TPM_BUCKET = pow(2, MAX_TPM_EXP);

    private static final double STRONG_BINDER_LIKELIHOOD = 0.001;

    public ExpressionDistribution(int totalPeptides)
    {
        mTpmBuckets = Lists.newArrayListWithCapacity(22);

        double tpmBucket = MIN_TPM_BUCKET;
        while(tpmBucket <= MAX_TPM_BUCKET)
        {
            mTpmBuckets.add(new TpmBucket(tpmBucket));
            tpmBucket *= 2;
        }

        mDecileBuckets = Lists.newArrayListWithCapacity(10);
        mRoundedTpmBuckets = Lists.newArrayList();
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

    public void formDistributions()
    {
        // form deciles
        int totalPeptides = mRoundedTpmBuckets.stream().mapToInt(x -> x.TotalCount).sum();
        int peptidesPerDecile = totalPeptides / 10;

        TpmBucket decileBucket = new TpmBucket(0.1);
        mDecileBuckets.add(decileBucket);

        for(int i = 0; i < mRoundedTpmBuckets.size(); ++i)
        {
            TpmBucket bucket = mRoundedTpmBuckets.get(i);

            if(decileBucket.TotalCount + bucket.TotalCount >= peptidesPerDecile && mDecileBuckets.size() < 10) //  || i == mAllExpData.size() - 1)
            {
                int required = peptidesPerDecile - decileBucket.TotalCount;
                int excess = bucket.TotalCount - required;
                double requiredFraction = required / (double)bucket.TotalCount;
                double excessFraction = 1 - requiredFraction;

                decileBucket.TotalCount += required;
                decileBucket.StrongBinders += (int)round(requiredFraction * bucket.StrongBinders);

                decileBucket = new TpmBucket(decileBucket.Bucket + 0.1);
                mDecileBuckets.add(decileBucket);

                decileBucket.TotalCount += excess;
                decileBucket.StrongBinders += (int)round(excessFraction * bucket.StrongBinders);
            }
            else
            {
                decileBucket.TotalCount += bucket.TotalCount;
                decileBucket.StrongBinders += bucket.StrongBinders;
            }
        }
    }

    private static double roundTpm(double tpm)
    {
        if(Doubles.equal(tpm, 0))
            return 0;

        // for deciles - round to 3 digits of precison
        int scale = (int)round(Math.log10(tpm));
        double discrete = pow(10, scale) / 100.0;

        return discrete * round(tpm / discrete);
    }

    private void addToRankedList(final PeptideExpressionData pepExpData)
    {
        // early exits
        double tpmRounded = roundTpm(pepExpData.tpm());

        int index = 0;

        while(index < mRoundedTpmBuckets.size())
        {
            TpmBucket bucket = mRoundedTpmBuckets.get(index);

            if(Doubles.equal(bucket.Bucket, tpmRounded))
            {
                bucket.add(pepExpData);
                return;
            }

            if(tpmRounded < bucket.Bucket)
                break;

            ++index;
        }

        TpmBucket bucket = new TpmBucket(tpmRounded);
        mRoundedTpmBuckets.add(index, bucket);
        bucket.add(pepExpData);
    }

    private void addToTpmBucket(final PeptideExpressionData pepExpData)
    {
        double rawTpm = pepExpData.tpm();
        double tpmRounded = pow(2, max(min(round(log(2, rawTpm)), MAX_TPM_EXP), MIN_TPM_EXP));  // 2**pmax(-5,pmin(14,round(log2(TPM),0));

        int index = 0;
        while(index < mTpmBuckets.size())
        {
            TpmBucket bucket = mTpmBuckets.get(index);

            if(Doubles.equal(bucket.Bucket, tpmRounded) || tpmRounded < bucket.Bucket)
                break;

            ++index;
        }

        TpmBucket bucket = mTpmBuckets.get(index);
        bucket.add(pepExpData);
    }

    public void writeDistributions(final String outputDir, final String outputId)
    {
        final String filename = formFilename(outputDir, FILE_ID_EXPRESSION_DIST, outputId);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Type,Bucket,BindingRate,TotalCount");
            writer.newLine();

            for(TpmBucket bucket : mTpmBuckets)
            {
                writer.write(String.format("%s,%4.3e,%.6f,%d",
                        EXP_TYPE_TPM_LEVEL, bucket.Bucket, bucket.calcRate(), bucket.TotalCount));

                writer.newLine();
            }

            for(TpmBucket bucket : mDecileBuckets)
            {
                writer.write(String.format("%s,%.1f,%.6f,%d",
                        EXP_TYPE_DECILE, bucket.Bucket, bucket.calcRate(), bucket.TotalCount));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write random peptide likelihood file: {}", e.toString());
        }
    }

    private class TpmBucket
    {
        public double Bucket;
        public int TotalCount;
        public int StrongBinders;

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

        public double calcRate()
        {
            int proteome = TotalCount - StrongBinders;
            return proteome > 0 ? StrongBinders / (double)proteome : 0;
        }

        public String toString() { return String.format("bucket(%.6f) total(%d) binders(%d)", Bucket, TotalCount, StrongBinders); }
    }
}
