package com.hartwig.hmftools.cup.somatics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sigs.PositionFrequencies.buildStandardChromosomeLengths;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.getChromosomeFromIndex;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.getPositionFromIndex;
import static com.hartwig.hmftools.common.sigs.SnvSigUtils.populateBucketMap;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_BUCKET_SIZE;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.extractPositionFrequencyCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_13;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_2;
import static com.hartwig.hmftools.cup.somatics.TrinucleotideCounts.extractTrinucleotideCounts;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.CategoryType;
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

    private final Map<String,Integer> mTriNucBucketNameMap;
    private final List<String> mSnv96BucketNames;
    private double[] mTriNucCounts;
    private int mTotalSnvCount = 0;

    private final SomaticSigs mSomaticSigs;

    private final PositionFrequencies mPosFrequencies;

    public SomaticVariantPrep(final PrepConfig config)
    {
        mConfig = config;

        mTriNucBucketNameMap = Maps.newHashMap();
        mSnv96BucketNames = Lists.newArrayList();
        populateBucketMap(mTriNucBucketNameMap, mSnv96BucketNames);

        mSomaticSigs = new SomaticSigs(null);

        // could add bucket size and max counts as config
        mPosFrequencies = new PositionFrequencies(
                mConfig.RefGenVersion,
                GEN_POS_BUCKET_SIZE,
                GEN_POS_MAX_SAMPLE_COUNT,
                buildStandardChromosomeLengths(mConfig.RefGenVersion),
                false
        );
    }

    @Override
    public CategoryType categoryType() { return CategoryType.SNV; }

    private void loadVariants(String sampleId)
    {
        List<SomaticVariant> variants = SomaticVariantsLoader.loadFromConfig(mConfig, sampleId, Lists.newArrayList(SNP));
        mVariants.addAll(variants);
    }

    private void getTrinucleotideCounts()
    {
        // build the 96 trinucleotide context counts
        mTriNucCounts = extractTrinucleotideCounts(mVariants, mTriNucBucketNameMap);

        for(int b = 0; b < mSnv96BucketNames.size(); ++b)
        {
            String bucketName = mSnv96BucketNames.get(b);
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
        // build genomic position counts
        AidApobecStatus aidApobecStatus = AidApobecStatus.FALSE_ONLY;
        extractPositionFrequencyCounts(mVariants, mPosFrequencies, aidApobecStatus);

        final int[] genPosCount = mPosFrequencies.getCounts();

        for(int b = 0; b < mPosFrequencies.getBucketCount(); ++b)
        {
            final String chromosome = getChromosomeFromIndex(mConfig.RefGenVersion, mPosFrequencies.chromosomePosIndex(), b);
            int position = getPositionFromIndex(mPosFrequencies.chromosomePosIndex(), chromosome, b, mPosFrequencies.getBucketSize());
            String keyName = format("%s_%d", chromosome, position);

            DataItem dataItem = new DataItem(DNA, ItemType.GEN_POS, keyName, genPosCount[b]);
            mDataItems.add(dataItem);
        }
    }

    private void getSignatureAllocations()
    {
        final double[] sigAllocations = mSomaticSigs.fitSampleCounts(mTriNucCounts);

        Map<String,Double> reportedAllocations = Maps.newHashMap();

        for(int i = 0; i < sigAllocations.length; ++i)
        {
            final String sigName = mSomaticSigs.getSigName(i);
            reportedAllocations.put(sigName, sigAllocations[i]);
        }

        for(Map.Entry<String,String> entry :SomaticSigs.REPORTABLE_SIGS.entrySet())
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
            CUP_LOGGER.error("sample({}) failed to extract somatic variant features: {}", sampleId, e.toString());
            System.exit(1);
        }

        return mDataItems;
    }
}
