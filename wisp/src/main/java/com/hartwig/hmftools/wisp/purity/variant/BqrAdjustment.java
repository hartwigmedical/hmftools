package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sage.SageCommon.generateBqrFilename;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.qual.BaseQualAdjustment;
import com.hartwig.hmftools.common.qual.BqrFile;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.qual.BqrRecord;
import com.hartwig.hmftools.wisp.purity.PurityConfig;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class BqrAdjustment
{
    private final PurityConfig mConfig;

    private final List<BqrContextData> mBqrContextData;

    public BqrAdjustment(final PurityConfig config)
    {
        mConfig = config;
        mBqrContextData = Lists.newArrayList();
    }

    private static final int BQR_MIN_QUAL = 35;

    private static final byte NO_KEY_VALUE = 1;

    private static final List<Integer> QUAL_THRESHOLDS = Lists.newArrayList(1, 38, 42, 45, 48);

    public void processSample(final String sampleId, final List<SomaticVariant> variants)
    {
        // calculate BQR error rates for each TNC and alt
        loadBqrData(sampleId);

        if(mBqrContextData.isEmpty())
            return;

        for(SomaticVariant variant : variants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            if(sampleFragData == null)
                continue;

            BqrContextData bqrContextData = getOrCreate(variant.decorator().trinucleotideContext(), variant.Alt);

            if(bqrContextData == null)
                continue;

            bqrContextData.SampleDepthTotal += sampleFragData.Depth;
            bqrContextData.SampleFragmentTotal += sampleFragData.AlleleCount;
        }

        // now test out each BQR qual threshold
        for(Integer qualThreshold : QUAL_THRESHOLDS)
        {
            checkQualThreshold(sampleId, qualThreshold);
        }
    }

    private void checkQualThreshold(final String sampleId, final int qualThreshold)
    {
        int depthTotal = 0;
        int fragmentTotal = 0;

        double expectedErrors = 0;

        for(BqrContextData bqrContextData : mBqrContextData)
        {
            if(bqrContextData.calculatedQual() < qualThreshold)
                continue;

            if(bqrContextData.SampleDepthTotal == 0)
                continue;

            depthTotal += bqrContextData.SampleDepthTotal;
            fragmentTotal += bqrContextData.SampleFragmentTotal;

            expectedErrors += bqrContextData.errorRatePerMillion() * bqrContextData.SampleDepthTotal / bqrContextData.SampleDepthTotal;

            // summarise(BQRErrorPerM=round(sum(BQRErrorPerM*SampleDP)/sum(SampleDP),0),SampleAD=sum(SampleAD),SampleDP=sum(SampleDP)) %>%
        }

        double errorRate = sampleErrorPerMillion(depthTotal, fragmentTotal);

        PoissonDistribution distribution = new PoissonDistribution(expectedErrors);

        double poisProbability  = distribution.cumulativeProbability((int)round(errorRate));

        CT_LOGGER.debug(format("sample(%s) qualThreshold(%d) frags(dp=%d ad=%d) errorRate(%.4f) prob(%.6f)",
                sampleId, qualThreshold, depthTotal, fragmentTotal, errorRate, poisProbability));

        /*

          summarise(BQRErrorPerM=round(sum(BQRErrorPerM*SampleDP)/sum(SampleDP),0),
            SampleAD=sum(SampleAD),
            SampleDP=sum(SampleDP)) %>%
  mutate(SampleErrorPerM=round(SampleAD/SampleDP*1e6,0),
         qualThreshold=40,
         pvalue=ppois(SampleAD,BQRErrorPerM*SampleDP/1e6,FALSE))

         */

    }

    public static double sampleErrorPerMillion(int depthTotal, int fragmentTotal)
    {
        return depthTotal > 0 ? round(100.0 * fragmentTotal / depthTotal * 1_000_000) / 100.0 : 0;
    }

    private void loadBqrData(final String sampleId)
    {
        String vcfDir = !mConfig.SomaticVcf.isEmpty() ? pathFromFile(mConfig.getSomaticVcf(sampleId)) : mConfig.SomaticDir;
        String bqrFilename = generateBqrFilename(vcfDir, sampleId);

        if(!Files.exists(Paths.get(bqrFilename)))
            return;

        List<BqrRecord> allCounts = BqrFile.read(bqrFilename);

        Map<BqrKey,Integer> summaryCounts = Maps.newHashMap();

        for(BqrRecord bqrRecord : allCounts)
        {
            BqrKey key = bqrRecord.Key;

            if(bqrRecord.Key.Quality <= BQR_MIN_QUAL)
                continue;

            if(bqrRecord.Key.ReadType == BqrReadType.DUAL)
                continue;

            BqrKey noAltKey = new BqrKey(key.Ref, NO_KEY_VALUE, key.TrinucleotideContext, key.Quality, key.ReadType);

            int count = summaryCounts.getOrDefault(noAltKey, 0);
            summaryCounts.put(noAltKey, count + bqrRecord.Count);

            if(key.Ref != key.Alt)
            {
                BqrContextData bqrErrorRate = getOrCreate(key.TrinucleotideContext, key.Alt);
                bqrErrorRate.AltCount += bqrRecord.Count;
            }
        }

        for(BqrContextData bqrContextData : mBqrContextData)
        {
            for(Map.Entry<BqrKey,Integer> entry : summaryCounts.entrySet())
            {
                BqrKey key = entry.getKey();

                if(new String(key.TrinucleotideContext).equals(bqrContextData.TrinucleotideContext))
                {
                    bqrContextData.TotalCount += entry.getValue();
                }
            }
        }

        for(BqrContextData bqrContextData : mBqrContextData)
        {
            CT_LOGGER.trace("sample({}) simple BQR summary: {}", sampleId, bqrContextData);
        }
    }

    private BqrContextData getOrCreate(final byte[] trinucleotideContext, final byte alt)
    {
        return getOrCreate(new String(trinucleotideContext), String.valueOf((char)alt));
    }

    private BqrContextData getOrCreate(final String trinucleotideContext, final String alt)
    {
        BqrContextData bqrErrorRate = mBqrContextData.stream()
                .filter(x -> x.Alt.equals(alt) && x.TrinucleotideContext.equals(trinucleotideContext))
                .findFirst().orElse(null);

        if(bqrErrorRate != null)
            return bqrErrorRate;

        bqrErrorRate = new BqrContextData(trinucleotideContext, alt);
        mBqrContextData.add(bqrErrorRate);
        return bqrErrorRate;
    }

    private class BqrContextData
    {
        public final String TrinucleotideContext;
        public final String Alt;

        public int AltCount;
        public int TotalCount;

        public int SampleFragmentTotal;
        public int SampleDepthTotal;

        public BqrContextData(final String trinucleotideContext, final String alt)
        {
            TrinucleotideContext = trinucleotideContext;
            Alt = alt;

            AltCount = 0;
            TotalCount = 0;

            SampleFragmentTotal = 0;
            SampleDepthTotal = 0;
        }

        private static final int PER_MILLION = 1_000_000;

        public double errorRatePerMillion()
        {
            return TotalCount > 0 ? round(100.0 * AltCount / TotalCount * 1_000_000) / 100.0 : 0;
        }

        public double calculatedQual()
        {
            return TotalCount > 0 ? BaseQualAdjustment.probabilityToPhredQual(AltCount / (double)TotalCount) : 0;
        }

        public String toString()
        {
            return format("key(%s:%s) counts(total=%d alt=%d) errorPerM(%.4f) calcQual(%.4f)",
                    TrinucleotideContext, Alt, TotalCount, AltCount, errorRatePerMillion(), calculatedQual());
        }
    }
}
