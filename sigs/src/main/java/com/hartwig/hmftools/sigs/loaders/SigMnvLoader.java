package com.hartwig.hmftools.sigs.loaders;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.sigs.common.CommonUtils.getNewFile;
import static com.hartwig.hmftools.sigs.loaders.SigSnvLoader.getBucketNameByIndex;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigMnvLoader
{
    private final DataLoaderConfig mConfig;

    private final Map<String,Integer> mBucketStringToIndex;
    private Matrix mSampleBucketCounts;

    private static final String MNV_3_BASES = "3Bases";
    private static final String MNV_4P_BASES = "4+Bases";

    private static final Logger LOGGER = LogManager.getLogger(SigMnvLoader.class);

    public SigMnvLoader(final DataLoaderConfig config)
    {
        mConfig = config;

        mBucketStringToIndex = Maps.newHashMap();
        buildBucketMap();
    }

    private void buildBucketMap()
    {
        mBucketStringToIndex.put(MNV_3_BASES, 0);
        mBucketStringToIndex.put(MNV_4P_BASES, 1);

        char[] bases = {'A','C', 'G', 'T'};

        List<String> basePairs = Lists.newArrayList();

        for(int i = 0; i < bases.length; ++i)
        {
            char base1 = bases[i];

            for (int j = 0; j < bases.length; ++j)
            {
                char base2 = bases[j];

                basePairs.add(String.format("%c%c", base1, base2));
            }
        }

        for(int i = 0; i < basePairs.size(); ++i)
        {
            final String bp1 = basePairs.get(i);

            for (int j = 0; j < basePairs.size(); ++j)
            {
                final String bp2 = basePairs.get(j);

                if(bp1.charAt(0) == bp2.charAt(0) || bp1.charAt(1) == bp2.charAt(1))
                    continue;

                String mutation = convertBasePair(bp1, bp2);

                if(mBucketStringToIndex.containsKey(mutation))
                    continue;

                mBucketStringToIndex.put(mutation, mBucketStringToIndex.size());
            }
        }
    }

    public void loadData(DatabaseAccess dbAccess)
    {
        mSampleBucketCounts = new Matrix(mBucketStringToIndex.size(), mConfig.SampleIds.size());

        LOGGER.info("retrieving MNV data for {} samples", mConfig.SampleIds.size());

        for(int sampleIndex = 0; sampleIndex < mConfig.SampleIds.size(); ++sampleIndex)
        {
            String sampleId = mConfig.SampleIds.get(sampleIndex);
            final List<SomaticVariant> variants = dbAccess.readSomaticVariants(sampleId, VariantType.MNP);

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

    private void processSampleVariants(final String sampleId, List<SomaticVariant> variants, int sampleIndex)
    {
        double[][] sampleCounts = mSampleBucketCounts.getData();

        for(final SomaticVariant variant : variants)
        {
            if(variant.isFiltered() || variant.type() != VariantType.MNP)
                continue;

            String rawContext = variant.trinucleotideContext();

            if(rawContext.contains("N"))
                continue;

            // check filters
            if(!mConfig.Filters.passesFilters(variant))
                continue;

            final String bucketName = mutationType(variant);

            Integer bucketIndex = mBucketStringToIndex.get(bucketName);

            if(bucketIndex == null)
            {
                LOGGER.error("sample({}) invalid bucketName({}) from ref({}) to alt({})",
                        sampleId, bucketName, variant.ref(), variant.alt());

                return;
            }

            ++sampleCounts[bucketIndex][sampleIndex];
        }
    }

    private String convertBasePair(final String ref, final String alt)
    {
        if(ref.equals("CA") || ref.equals("AA") || ref.equals("AG") || ref.equals("GA") || ref.equals("GG") || ref.equals("GT"))
        {
            return String.format("%c%c>%c%c",
                    swapDnaBase(ref.charAt(1)), swapDnaBase(ref.charAt(0)), swapDnaBase(alt.charAt(1)), swapDnaBase(alt.charAt(0)));

            // paste(baseConvert(substr(ref,2,2)), baseConvert(substr(ref,1,1)), '>', baseConvert(substr(alt,2,2)), baseConvert(substr(alt,1,1)), sep='')
        }

        if(ref.equals("AT") || ref.equals("TA") || ref.equals("CG") || ref.equals("GC"))
        {
            // these refs will all remain the same with a reversal, so then convert the alts
            if(alt.equals("CA") || alt.equals("AA") || alt.equals("AG") || alt.equals("GA") || alt.equals("GG") || alt.equals("GT"))
            {
                // mutation = paste(ref, '>', baseConvert(substr(alt,2,2)), baseConvert(substr(alt,1,1)), sep='')
                return String.format("%s>%c%c", ref, swapDnaBase(alt.charAt(1)), swapDnaBase(alt.charAt(0)));
            }
        }

        return ref + ">" + alt;
    }

    private String mutationType(final SomaticVariant variant)
    {
        int altLength = variant.alt().length();

        if(altLength == 2)
        {
            return convertBasePair(variant.ref(), variant.alt());
        }
        else if(altLength == 3)
        {
            return MNV_3_BASES;
        }
        else
        {
            return MNV_4P_BASES;
        }
    }

}
