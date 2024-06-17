package com.hartwig.hmftools.sigs.loaders;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.sigs.SnvSigUtils.contextFromVariant;
import static com.hartwig.hmftools.common.sigs.SnvSigUtils.populateBucketMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class SigSnvLoader
{
    private final VariantFilters mFilters;
    private final List<String> mSampleIds;

    private final Map<String,Integer> mBucketStringToIndex;
    private Matrix mSampleBucketCounts;

    private final List<PositionFrequencies> mPositionFrequencies;

    private final List<BufferedWriter> mPosFreqWriters;

    private static final int SNV_BUCKET_COUNT = 96;

    public SigSnvLoader(final VariantFilters filters)
    {
        mFilters = filters;
        mSampleIds = Lists.newArrayList();

        mBucketStringToIndex = Maps.newHashMap();
        populateBucketMap(mBucketStringToIndex);

        mPositionFrequencies = Lists.newArrayList();
        mPosFreqWriters = Lists.newArrayList();
    }

    public void setSampleIds(final List<String> sampleIds)
    {
        mSampleIds.clear();
        mSampleIds.addAll(sampleIds);
    }

    public void initialisePositionFrequencies(final String outputDir, final List<Integer> bucketSizes)
    {
        mPositionFrequencies.clear();

        bucketSizes.forEach(x -> mPositionFrequencies.add(new PositionFrequencies(V37, x)));

        mPositionFrequencies.forEach(x -> mPosFreqWriters.add(PositionFrequencies.createFrequencyCountsWriter(outputDir, x.getBucketSize())));
    }

    public final List<PositionFrequencies> getPositionFrequencies() { return mPositionFrequencies; }

    public Matrix getSampleBucketCounts() { return mSampleBucketCounts; }

    public void loadData(final DatabaseAccess dbAccess, final String vcfFile, boolean writePosFreqData)
    {
        mSampleBucketCounts = new Matrix(SNV_BUCKET_COUNT, mSampleIds.size());

        if(mSampleIds.size() > 1)
        {
            SIG_LOGGER.debug("retrieving SNV data for {} samples", mSampleIds.size());
        }

        for(int sampleIndex = 0; sampleIndex < mSampleIds.size(); ++sampleIndex)
        {
            String sampleId = mSampleIds.get(sampleIndex);

            final List<SomaticVariant> variants = vcfFile != null ?
                    loadSomaticVariants(vcfFile, sampleId) : dbAccess.readSomaticVariants(sampleId, VariantType.SNP);

            SIG_LOGGER.info("sample({}) processing {} variants", sampleId, variants.size());

            processSampleVariants(sampleId, variants, sampleIndex);

            if(writePosFreqData)
            {
                for(int i = 0; i < mPositionFrequencies.size(); ++i)
                {
                    final PositionFrequencies positionFrequency = mPositionFrequencies.get(i);
                    positionFrequency.writeFrequencyCounts(mPosFreqWriters.get(i), sampleId);
                    positionFrequency.clear();
                }
            }

            if(sampleIndex > 0 && (sampleIndex % 100) == 0)
            {
                SIG_LOGGER.info("processed {} samples", sampleIndex);
            }
        }

        if(writePosFreqData)
            mPosFreqWriters.forEach(x -> closeBufferedWriter(x));
    }

    private List<SomaticVariant> loadSomaticVariants(final String vcfFile, final String sampleId)
    {
        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        SomaticVariantFactory variantFactory = new SomaticVariantFactory(filter);
        final List<SomaticVariant> variantList = Lists.newArrayList();

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        for(VariantContext variant : vcfFileReader.iterator())
        {
            if(variant.isFiltered())
                continue;

            final SomaticVariant somaticVariant = variantFactory.createVariant(sampleId, variant).orElse(null);

            if(somaticVariant == null || !somaticVariant.isSnp())
                continue;

            variantList.add(somaticVariant);
        }

        return variantList;
    }

    public void writeSampleCounts(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("BucketName");

            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                writer.write(String.format(",%s", mSampleIds.get(i)));
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
            SIG_LOGGER.error("error writing to outputFile: {}", e.toString());
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
            if(mFilters != null && !mFilters.passesFilters(variant))
                continue;

            for(PositionFrequencies positionFrequencies : mPositionFrequencies)
            {
                positionFrequencies.addPosition(variant.chromosome(), (int)variant.position());
            }

            final String bucketName = contextFromVariant(variant);
            Integer bucketIndex = mBucketStringToIndex.get(bucketName);

            if(bucketIndex == null)
            {
                SIG_LOGGER.error("sample({}) invalid bucketName({}) from var({}>{}) context={})",
                        sampleId, bucketName, variant.ref(), variant.alt(), variant.trinucleotideContext());

                return;
            }

            ++sampleCounts[bucketIndex][sampleIndex];
        }
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

}
