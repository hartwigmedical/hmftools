package com.hartwig.hmftools.cup.traits;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.GENDER;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.MS_INDELS_TMB;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.WGD;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class SampleTraitPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    private static final String FLOAT_FORMAT_MS_INDELS_TMB = "%.4f";

    public SampleTraitPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return CategoryType.SAMPLE_TRAIT; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        try
        {
            final PurityContext purityContext = PurityContextFile.readWithQC(
                    mConfig.purpleQcFile(sampleId),
                    mConfig.purplePurityFile(sampleId));

            dataItems.add(new DataItem(
                    DNA, ItemType.SAMPLE_TRAIT, GENDER.getAlias(),
                    purityContext.gender() == Gender.MALE));

            dataItems.add(new DataItem(
                    DNA, ItemType.TUMOR_MUTATIONAL_BURDEN, MS_INDELS_TMB.getAlias(),
                    purityContext.microsatelliteIndelsPerMb(), FLOAT_FORMAT_MS_INDELS_TMB));

            dataItems.add(new DataItem(
                    DNA, ItemType.SAMPLE_TRAIT, WGD.getAlias(),
                    purityContext.wholeGenomeDuplication()));

            return dataItems;
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to extract category({}):", sampleId, categoryType());
            e.printStackTrace();
            System.exit(1);
        }

        return dataItems;
    }
}
