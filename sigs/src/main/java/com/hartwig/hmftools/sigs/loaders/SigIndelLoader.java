package com.hartwig.hmftools.sigs.loaders;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.sigs.fitter.VariantFilters;

public class SigIndelLoader
{
    private final VariantFilters mFilters;

    private final Map<String,Integer> mBucketStringToIndex;
    private Matrix mSampleBucketCounts;

    private static final int MAX_GROUP_COUNT = 4;

    public SigIndelLoader(final VariantFilters filters)
    {
        mFilters = filters;

        mBucketStringToIndex = Maps.newHashMap();
        buildBucketMap();
    }

    private void buildBucketMap()
    {
        int index = 0;

        for(int i = 1; i <= MAX_GROUP_COUNT; ++i)
        {
            mBucketStringToIndex.put(String.format("DEL_%d_MH_RCH", i), index++);
            mBucketStringToIndex.put(String.format("DEL_%d_NoMH_RCL", i), index++);
            mBucketStringToIndex.put(String.format("DEL_%d_MH_RCL", i), index++);
        }

        for(int i = 1; i <= MAX_GROUP_COUNT; ++i)
        {
            mBucketStringToIndex.put(String.format("INS_%d_NoMH_RCL", i), index++);
            mBucketStringToIndex.put(String.format("INS_%d_NoMH_RCH", i), index++);
        }
    }

    private void processSampleVariants(final String sampleId, List<SmallVariant> variants, int sampleIndex)
    {
        double[][] sampleCounts = mSampleBucketCounts.getData();

        for(final SmallVariant variant : variants)
        {
            if(variant.isFiltered())
                continue;

            // check filters
            if(!mFilters.passesFilters(variant))
                continue;

            final String bucketName = convertToBucket(variant);
            Integer bucketIndex = mBucketStringToIndex.get(bucketName);

            if(bucketIndex == null)
            {
                SIG_LOGGER.error("sample({}) invalid bucketName({}) from ref({}) alt({})",
                        sampleId, bucketName, variant.ref(), variant.alt());

                return;
            }

            ++sampleCounts[bucketIndex][sampleIndex];
        }
    }

    private static final String SUBTYPE_DEL = "DEL";
    private static final String SUBTYPE_INS = "INS";

    private String convertToBucket(final SmallVariant variant)
    {
        int refLength = variant.ref().length();
        int altLength = variant.alt().length();

        if(variant.alt().contains(","))
        {
            // take the first alt's length if there are more than 1
            String[] altSplits = variant.alt().split(",");
            altLength = altSplits[0].length();
        }

        final String subtype = altLength > refLength ? SUBTYPE_INS : SUBTYPE_DEL;

        int length = abs(altLength - refLength);
        int lengthGroup = min(length, MAX_GROUP_COUNT);
        boolean highRepeatCount = variant.repeatCount() >= 4;

        boolean hasMicrohomology = !variant.microhomology().isEmpty() && !variant.microhomology().equals(".");

        boolean microhomology = subtype == SUBTYPE_DEL && highRepeatCount && hasMicrohomology;

        // make a combined bucket name
        return String.format("%s_%d_%s_%s", subtype, lengthGroup, microhomology ? "MH" : "NoMH", highRepeatCount ? "RCH" : "RCL");
    }

}
