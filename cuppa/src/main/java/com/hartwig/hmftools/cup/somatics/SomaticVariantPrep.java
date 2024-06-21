package com.hartwig.hmftools.cup.somatics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_BUCKET_SIZE;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_13;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_2;

import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SnvSigUtils;
import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;
import com.hartwig.hmftools.cup.traits.SampleTraitType;

public class SomaticVariantPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    private static final String FLOAT_FORMAT_SIG_ALLOCATION = "%.1f";

    private final List<SomaticVariant> mVariants = new ArrayList<>();
    private final List<DataItem> mDataItems = new ArrayList<>();

    private double[] mTriNucCounts;
    private int mTotalSnvCount = 0;

    public SomaticVariantPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return CategoryType.SNV; }

    private void loadVariants(String sampleId) throws NoSuchFileException
    {
        List<SomaticVariant> variants = SomaticVariantsLoader.loadFromConfig(mConfig, sampleId, Lists.newArrayList(SNP));
        mVariants.addAll(variants);
    }

    private void getTrinucleotideCounts()
    {
        // build the 96 trinucleotide context counts
        Map<String,Integer> triNucBucketNameMap = new HashMap<>();
        List<String> snv96BucketNames = new ArrayList<>();
        SnvSigUtils.populateBucketMap(triNucBucketNameMap, snv96BucketNames);
        mTriNucCounts = TrinucleotideCounts.extractTrinucleotideCounts(mVariants, triNucBucketNameMap);

        for(int b = 0; b < snv96BucketNames.size(); ++b)
        {
            String bucketName = snv96BucketNames.get(b);
            int count = (int) mTriNucCounts[b];

            mTotalSnvCount += count;

            DataItem dataItem = new DataItem(DNA, ItemType.SNV96, bucketName, count);
            mDataItems.add(dataItem);
        }
    }

    private void getSnvCount()
    {
        DataItem dataItem = new DataItem(DNA, ItemType.TUMOR_MUTATIONAL_BURDEN, SampleTraitType.SNV_COUNT.getAlias(), mTotalSnvCount);
        mDataItems.add(dataItem);
    }

    private void getGenomicPositionCounts()
    {
        PositionFrequencies posFrequencies = new PositionFrequencies(
                mConfig.RefGenVersion,
                GEN_POS_BUCKET_SIZE,
                PositionFrequencies.buildStandardChromosomeLengths(mConfig.RefGenVersion),
                false
        );

        // build genomic position counts
        AidApobecStatus aidApobecStatus = AidApobecStatus.FALSE_ONLY;
        GenomicPositions.extractPositionFrequencyCounts(mVariants, posFrequencies, aidApobecStatus);

        final int[] genPosCount = posFrequencies.getCounts();

        String chromosomeY = mConfig.RefGenVersion.versionedChromosome("chrY");

        for(int b = 0; b < posFrequencies.getBucketCount(); ++b)
        {
            final String chromosome = PositionFrequencies.getChromosomeFromIndex(mConfig.RefGenVersion, posFrequencies.chromosomePosIndex(), b);

            if(chromosome.equals(chromosomeY))
            {
                // gen_pos frequencies (i.e. SNV counts) of chromosome Y correlate directly with sex, but this is already covered by the
                // 'trait.is_male' feature. Ignore chromosome Y to prevent duplicating features.
                continue;
            }

            int position = PositionFrequencies.getPositionFromIndex(posFrequencies.chromosomePosIndex(), chromosome, b, posFrequencies.getBucketSize());
            String keyName = format("%s_%d", chromosome, position);

            DataItem dataItem = new DataItem(DNA, ItemType.GEN_POS, keyName, genPosCount[b]);
            mDataItems.add(dataItem);
        }
    }

    private void getSignatureAllocations()
    {
        SomaticSigs somaticSigs = new SomaticSigs(null);

        final double[] sigAllocations = somaticSigs.fitSampleCounts(mTriNucCounts);

        Map<String,Double> reportedAllocations = Maps.newHashMap();

        for(int i = 0; i < sigAllocations.length; ++i)
        {
            final String sigName = somaticSigs.getSigName(i);
            reportedAllocations.put(sigName, sigAllocations[i]);
        }

        for(Map.Entry<String,String> entry : SomaticSigs.REPORTABLE_SIGS.entrySet())
        {
            String sigName = entry.getKey();
            double sigAllocation = reportedAllocations.get(sigName);

            // combine 2 & 13
            if(sigName.equalsIgnoreCase(SIG_NAME_13))
                continue;

            if(sigName.equalsIgnoreCase(SIG_NAME_2))
            {
                sigAllocation += reportedAllocations.get(SIG_NAME_13);
            }

            DataItem dataItem = new DataItem(DNA, ItemType.SIGNATURE, entry.getValue(), sigAllocation, FLOAT_FORMAT_SIG_ALLOCATION);
            mDataItems.add(dataItem);
        }
    }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        try
        {
            loadVariants(sampleId);
            getTrinucleotideCounts();
            getSnvCount();
            getSignatureAllocations();
            getGenomicPositionCounts();
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to extract category({}):", sampleId, categoryType());
            e.printStackTrace();
            System.exit(1);
        }

        return mDataItems;
    }
}
