package com.hartwig.hmftools.sig_analyser.loaders;

import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.getNewFile;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.common.sigs.SigMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigSnvLoader
{
    private final DataLoaderConfig mConfig;

    private final Map<String,Integer> mBucketStringToIndex;
    private SigMatrix mSampleBucketCounts;

    private static final int SNV_BUCKET_COUNT = 96;

    private static final Logger LOGGER = LogManager.getLogger(SigSnvLoader.class);

    public SigSnvLoader(final DataLoaderConfig config)
    {
        mConfig = config;

        mBucketStringToIndex = Maps.newHashMap();
        buildBucketMap();
    }

    private void buildBucketMap()
    {
        char[] refBases = {'C', 'T'};
        char[] bases = {'A','C', 'G', 'T'};
        int index = 0;

        for(int i = 0; i < refBases.length; ++i)
        {
            char ref = refBases[i];

            for(int j = 0; j < bases.length; ++j)
            {
                char alt = bases[j];

                if(ref != alt)
                {
                    String baseChange = String.format("%c>%c", ref, alt);

                    for (int k = 0; k < bases.length; ++k)
                    {
                        char before = bases[k];

                        for (int l = 0; l < bases.length; ++l)
                        {
                            char after = bases[l];

                            String context = String.format("%c%c%c", before, ref, after);

                            String bucketName = baseChange + "_" + context;

                            mBucketStringToIndex.put(bucketName, index);
                            ++index;
                        }
                    }
                }
            }
        }
    }

    public void loadData(DatabaseAccess dbAccess)
    {
        mSampleBucketCounts = new SigMatrix(SNV_BUCKET_COUNT, mConfig.SampleIds.size());

        LOGGER.info("retrieving SNV data for {} samples", mConfig.SampleIds.size());

        for(int sampleIndex = 0; sampleIndex < mConfig.SampleIds.size(); ++sampleIndex)
        {
            String sampleId = mConfig.SampleIds.get(sampleIndex);
            final List<SomaticVariant> variants = dbAccess.readSomaticVariants(sampleId, VariantType.SNP);

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

            /*
            int i = 0;
            for(; i < mConfig.SampleIds.size()-1; ++i)
            {
                writer.write(String.format("%s,", mConfig.SampleIds.get(i)));
            }
            writer.write(String.format("%s", mConfig.SampleIds.get(i)));

            writer.newLine();

            writeMatrixData(writer, mSampleBucketCounts, true);
            */

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
            if(variant.isFiltered() || !variant.isSnp())
                continue;

            if(variant.alt().length() != 1)
                continue;

            String rawContext = variant.trinucleotideContext();

            if(rawContext.contains("N"))
                continue;

            // check filters
            if(!mConfig.passesFilters(variant))
                continue;

            // convert base change to standard set and the context accordingly
            String baseChange;
            String context;
            if(variant.ref().charAt(0) == 'A' || variant.ref().charAt(0) == 'G')
            {
                baseChange = String.format("%c>%c", convertBase(variant.ref().charAt(0)), convertBase(variant.alt().charAt(0)));

                // convert the context as well
                context = String.format("%c%c%c",
                        convertBase(rawContext.charAt(2)), convertBase(rawContext.charAt(1)), convertBase(rawContext.charAt(0)));
            }
            else
            {
                baseChange = variant.ref() + ">" + variant.alt();
                context = rawContext;
            }

            String bucketName = baseChange + "_" + context;
            Integer bucketIndex = mBucketStringToIndex.get(bucketName);

            if(bucketIndex == null)
            {
                LOGGER.error("sample({}) invalid bucketName({}) from baseChange({} raw={}>{}) context({} raw={}",
                        sampleId, bucketName, baseChange, variant.ref(), variant.alt(), context, rawContext);

                return;
            }

            ++sampleCounts[bucketIndex][sampleIndex];
        }
    }

    public static char convertBase(char base)
    {
        if(base == 'A') return 'T';
        if(base == 'T') return 'A';
        if(base == 'C') return 'G';
        if(base == 'G') return 'C';
        return base;
    }

    public static String getBucketNameByIndex(final Map<String,Integer> bucketNameIndexMap, int index)
    {
        for(Map.Entry<String,Integer> entry : bucketNameIndexMap.entrySet())
        {
            if(entry.getValue() == index)
                return entry.getKey();
        }

        return String.format("MissingBucket_%d", index);
    }

    private static String standardiseSnv(final String snv)
    {
        // convert to equivalent strand's base
        if(snv.equals("G>T")) return "C>A";
        if(snv.equals("G>C")) return "C>G";
        if(snv.equals("G>A")) return "C>T";
        if(snv.equals("A>T")) return "T>A";
        if(snv.equals("A>G")) return "T>C";
        if(snv.equals("A>C")) return "T>G";

        return snv;
    }

}
