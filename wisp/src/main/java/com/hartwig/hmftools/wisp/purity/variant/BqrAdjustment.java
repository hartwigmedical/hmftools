package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sage.SageCommon.generateBqrFilename;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.BQR_MIN_ERROR_RATE;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
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

    private static final byte NO_KEY_VALUE = 1;

    public boolean hasValidData() { return !mBqrContextData.isEmpty(); }

    public List<BqrContextData> getThresholdBqrData(final int qualThreshold)
    {
        return qualThreshold <= 0 ?
                mBqrContextData : mBqrContextData.stream().filter(x -> x.calculatedQual() >= qualThreshold).collect(Collectors.toList());
    }

    public double calcErrorRate(final String triNucContext, final String alt)
    {
        BqrContextData bqrData = mBqrContextData.stream()
                .filter(x -> x.TrinucleotideContext.equals(triNucContext) && x.Alt.equals(alt)).findFirst().orElse(null);

        return bqrData != null ? max(bqrData.errorRate(), BQR_MIN_ERROR_RATE) : BQR_MIN_ERROR_RATE;
    }

    public static double calcErrorRate(final List<BqrContextData> bqrContextData)
    {
        int depthTotal = 0;
        int fragmentTotal = 0;

        Set<String> processedTriNucContext = Sets.newHashSet();

        for(BqrContextData bqrData : bqrContextData)
        {
            if(bqrData.TotalCount == 0)
                continue;

            fragmentTotal += bqrData.AltCount;

            if(!processedTriNucContext.contains(bqrData.TrinucleotideContext))
            {
                processedTriNucContext.add(bqrData.TrinucleotideContext);
                depthTotal += bqrData.TotalCount;
            }
        }

        return depthTotal > 0 ? fragmentTotal / (double)depthTotal : 0;
    }

    public static boolean hasVariantContext(
            final List<BqrContextData> bqrContextData, final String trinucleotideContext, final String alt)
    {
        return bqrContextData.stream().anyMatch(x -> x.Alt.equals(alt) && x.TrinucleotideContext.equals(trinucleotideContext));
    }

    public void loadBqrData(final String sampleId)
    {
        mBqrContextData.clear();

        String bqrFileDir;

        if(mConfig.BqrDir != null)
            bqrFileDir = mConfig.BqrDir;
        else if(!mConfig.SomaticVcf.isEmpty())
            bqrFileDir = pathFromFile(mConfig.getSomaticVcf(sampleId));
        else
            bqrFileDir = mConfig.SomaticDir;

        String bqrFilename = generateBqrFilename(bqrFileDir, sampleId);

        if(!Files.exists(Paths.get(bqrFilename)))
        {
            CT_LOGGER.warn("sample({}) missing BQR file: {}", sampleId, bqrFilename);
            return;
        }

        List<BqrRecord> allCounts = BqrFile.read(bqrFilename);

        Map<BqrKey,Integer> summaryCounts = Maps.newHashMap();

        for(BqrRecord bqrRecord : allCounts)
        {
            BqrKey key = bqrRecord.Key;

            if(bqrRecord.Key.Quality < mConfig.BqrQualThreshold)
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
}
