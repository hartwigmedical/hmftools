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
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.qual.BqrFile;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.qual.BqrRecord;
import com.hartwig.hmftools.wisp.purity.PurityConfig;

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

    public List<BqrContextData> getThresholdBqrData(final double qualThreshold)
    {
        return qualThreshold <= 0 ?
                mBqrContextData : mBqrContextData.stream().filter(x -> x.calculatedQual() >= qualThreshold).collect(Collectors.toList());
    }

    public static double calcErrorRate(final List<BqrContextData> bqrContextData)
    {
        int depthTotal = 0;
        int fragmentTotal = 0;

        for(BqrContextData bqrData : bqrContextData)
        {
            if(bqrData.TotalCount == 0)
                continue;

            depthTotal += bqrData.TotalCount;
            fragmentTotal += bqrData.AltCount;
        }

        return depthTotal > 0 ? fragmentTotal / (double)depthTotal : 0;
        // return sampleErrorPerMillion(depthTotal, fragmentTotal);
    }

    public static boolean hasVariantContext(
            final List<BqrContextData> bqrContextData, final String trinucleotideContext, final String alt)
    {
        return bqrContextData.stream().anyMatch(x -> x.Alt.equals(alt) && x.TrinucleotideContext.equals(trinucleotideContext));
    }

    public static double sampleErrorPerMillion(int depthTotal, int fragmentTotal)
    {
        return depthTotal > 0 ? round(100.0 * fragmentTotal / depthTotal * 1_000_000) / 100.0 : 0;
    }

    public void loadBqrData(final String sampleId)
    {
        mBqrContextData.clear();

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

    /*
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
    */
}
