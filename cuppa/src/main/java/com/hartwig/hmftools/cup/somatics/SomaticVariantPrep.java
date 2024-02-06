package com.hartwig.hmftools.cup.somatics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sigs.PositionFrequencies.buildStandardChromosomeLengths;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.getChromosomeFromIndex;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.getPositionFromIndex;
import static com.hartwig.hmftools.common.sigs.SnvSigUtils.populateBucketMap;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.DATA_TYPE_SNV_COUNT;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_BUCKET_SIZE;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.extractPositionFrequencyCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_13;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_2;
import static com.hartwig.hmftools.cup.somatics.TrinucleotideCounts.extractTrinucleotideCounts;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class SomaticVariantPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    private final Map<String,Integer> mTriNucBucketNameMap;
    private final List<String> mSnv96BucketNames;

    private final PositionFrequencies mPosFrequencies;

    private final SomaticSigs mSomaticSigs;

    public SomaticVariantPrep(final PrepConfig config)
    {
        mConfig = config;

        mTriNucBucketNameMap = Maps.newHashMap();
        mSnv96BucketNames = Lists.newArrayList();
        populateBucketMap(mTriNucBucketNameMap, mSnv96BucketNames);

        mSomaticSigs = new SomaticSigs(null);

        // could add bucket size and max counts as config
        mPosFrequencies = new PositionFrequencies(
                mConfig.RefGenVersion, GEN_POS_BUCKET_SIZE, GEN_POS_MAX_SAMPLE_COUNT,
                buildStandardChromosomeLengths(mConfig.RefGenVersion), false);
    }

    @Override
    public CategoryType categoryType() { return CategoryType.SNV; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        final String purpleDataDir = mConfig.getPurpleDataDir(sampleId);

        try
        {
            List<SomaticVariant> variants = Lists.newArrayList();

            // load variants
            final String somaticVcfFile = PurpleCommon.purpleSomaticVcfFile(purpleDataDir, sampleId);
            variants.addAll(SomaticDataLoader.loadSomaticVariantsFromVcf(somaticVcfFile, Lists.newArrayList(SNP)));

            // build the 96 trinucleotide context counts
            final double[] triNucCounts = extractTrinucleotideCounts(variants, mTriNucBucketNameMap);

            int totalSnvCount = 0;

            for(int b = 0; b < mSnv96BucketNames.size(); ++b)
            {
                String bucketName = mSnv96BucketNames.get(b);
                int count = (int)triNucCounts[b];

                totalSnvCount += count;

                dataItems.add(new DataItem(DNA, ItemType.SNV96, bucketName, String.valueOf(count)));
            }

            dataItems.add(new DataItem(DNA, ItemType.SAMPLE_TRAIT, DATA_TYPE_SNV_COUNT, String.valueOf(totalSnvCount)));

            // build genomic position counts
            AidApobecStatus aidApobecStatus = AidApobecStatus.FALSE_ONLY;
            extractPositionFrequencyCounts(variants, mPosFrequencies, aidApobecStatus);
            final int[] genPosCount = mPosFrequencies.getCounts();

            for(int b = 0; b < mPosFrequencies.getBucketCount(); ++b)
            {
                final String chromosome = getChromosomeFromIndex(mConfig.RefGenVersion, mPosFrequencies.chromosomePosIndex(), b);
                int position = getPositionFromIndex(mPosFrequencies.chromosomePosIndex(), chromosome, b, mPosFrequencies.getBucketSize());
                String keyName = format("%s_%d", chromosome, position);

                dataItems.add(new DataItem(DNA, ItemType.GEN_POS, keyName, String.valueOf(genPosCount[b])));
            }

            // write signatures
            final double[] sigAllocations = mSomaticSigs.fitSampleCounts(triNucCounts);

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

                dataItems.add(new DataItem(DNA, ItemType.SIGNATURE, entry.getValue(), format("%.1f", sigAllocation)));
            }

            return dataItems;
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) sample traits - failed to load purity file from dir{}): {}",
                    sampleId, purpleDataDir, e.toString());

            return null;
        }
    }
}
