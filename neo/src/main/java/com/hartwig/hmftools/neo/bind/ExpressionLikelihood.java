package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.EXP_TYPE_TPM_LEVEL;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_DATA_TYPE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_TPM_BUCKET;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_TPM_RATE;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

public class ExpressionLikelihood
{
    private final List<TpmRate> mTpmRates;

    public static final String EXP_LIKELIHOOD_FILE = "exp_likelihood_file";

    public ExpressionLikelihood()
    {
        mTpmRates = Lists.newArrayList();
    }

    public boolean hasData() { return !mTpmRates.isEmpty(); }

    public double calcLikelihood(double tpm)
    {
        if(tpm < mTpmRates.get(0).TpmBucket)
            return mTpmRates.get(0).Likelihood;

        if(tpm > mTpmRates.get(mTpmRates.size() - 1).TpmBucket)
            return mTpmRates.get(mTpmRates.size() - 1).Likelihood;

        for(int i = 0; i < mTpmRates.size(); ++i)
        {
            TpmRate bucket = mTpmRates.get(i);

            if(Doubles.equal(tpm, bucket.TpmBucket))
                return bucket.Likelihood;

            TpmRate nextBucket = i < mTpmRates.size() - 1 ? mTpmRates.get(i + 1) : null;

            if(nextBucket != null && Doubles.equal(tpm, nextBucket.TpmBucket))
                return nextBucket.Likelihood;

            if(tpm > bucket.TpmBucket)
            {
                if(nextBucket == null)
                    break;

                if(tpm < nextBucket.TpmBucket)
                {
                    // interpolate between the distribution to set the rank
                    double upperPerc = (tpm - bucket.TpmBucket) / (nextBucket.TpmBucket - bucket.TpmBucket);
                    return upperPerc * nextBucket.Likelihood + (1 - upperPerc) * bucket.Likelihood;
                }
            }
        }

        return INVALID_CALC;
    }

    public boolean loadTpmRates(final String filename)
    {
        if(filename == null)
            return true;

        try
        {
            // Type,Bucket,BindingRate,TotalCount
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            int dataTypeIndex = fieldsIndexMap.get(FLD_DATA_TYPE);
            int bucketIndex = fieldsIndexMap.get(FLD_TPM_BUCKET);
            int rateIndex = fieldsIndexMap.get(FLD_TPM_RATE);
            
            for(String line : lines)
            {
                String[] items = line.split(DELIM, -1);

                if(!items[dataTypeIndex].equals(EXP_TYPE_TPM_LEVEL))
                    continue;

                double tpmBucket = Double.parseDouble(items[bucketIndex]);
                double tpmRate = Double.parseDouble(items[rateIndex]);
                
                mTpmRates.add(new TpmRate(tpmBucket, tpmRate));
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load expression likelhood rates file: {}" ,e.toString());
            return false;
        }

        return true;
    }

    private class TpmRate
    {
        public final double TpmBucket;
        public final double Likelihood;

        public TpmRate(final double tpmBucket, final double likelihood)
        {
            TpmBucket = tpmBucket;
            Likelihood = likelihood;
        }

        public String toString() { return String.format("TPM bucket(%.4f) likelihood(%.4f)", TpmBucket, Likelihood); }
    }
}
