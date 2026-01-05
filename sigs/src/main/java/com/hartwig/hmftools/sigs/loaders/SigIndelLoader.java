package com.hartwig.hmftools.sigs.loaders;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sigs.common.CommonUtils.getNewFile;
import static com.hartwig.hmftools.sigs.loaders.SigSnvLoader.getBucketNameByIndex;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigIndelLoader
{
    private final DataLoaderConfig mConfig;

    private final Map<String,Integer> mBucketStringToIndex;
    private Matrix mSampleBucketCounts;

    private static final int MAX_GROUP_COUNT = 4;

    private static final Logger LOGGER = LogManager.getLogger(SigIndelLoader.class);

    public SigIndelLoader(final DataLoaderConfig config)
    {
        mConfig = config;

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

    public void loadData(DatabaseAccess dbAccess)
    {
        mSampleBucketCounts = new Matrix(mBucketStringToIndex.size(), mConfig.SampleIds.size());

        LOGGER.info("retrieving INDEL data for {} samples", mConfig.SampleIds.size());

        for(int sampleIndex = 0; sampleIndex < mConfig.SampleIds.size(); ++sampleIndex)
        {
            String sampleId = mConfig.SampleIds.get(sampleIndex);
            final List<SmallVariant> variants = dbAccess.readSomaticVariants(sampleId, VariantType.INDEL);

            LOGGER.info("sample({}:{}) processing {} variants", sampleIndex, sampleId, variants.size());

            processSampleVariants(sampleId, variants, sampleIndex);
        }

        try
        {
            BufferedWriter writer = getNewFile(mConfig.OutputDir, mConfig.OutputFileId + "_sample_counts.csv");

            writer.write("BucketName");

            for(int i = 0; i < mConfig.SampleIds.size(); ++i)
            {
                writer.write(String.format(",%s", mConfig.SampleIds.get(i)));
            }

            writer.newLine();

            double[][] scData = mSampleBucketCounts.getData();

            for(int i = 0; i < mSampleBucketCounts.Rows; ++i)
            {
                writer.write(getBucketNameByIndex(mBucketStringToIndex, i));

                for(int j = 0; j < mSampleBucketCounts.Cols; ++j) {

                    writer.write(String.format(",%.0f", scData[i][j]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
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
            if(!mConfig.Filters.passesFilters(variant))
                continue;

            final String bucketName = convertToBucket(variant);
            Integer bucketIndex = mBucketStringToIndex.get(bucketName);

            if(bucketIndex == null)
            {
                LOGGER.error("sample({}) invalid bucketName({}) from ref({}) alt({})",
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
