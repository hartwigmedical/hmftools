package com.hartwig.hmftools.cup.traits;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.GENDER;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.MS_INDELS_TMB;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.WGD;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
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

    public SampleTraitPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return CategoryType.SAMPLE_TRAIT; }

    private static String boolToIntString(final boolean boolValue) { return boolValue ? "1" : "0"; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        try
        {
            final PurityContext purityContext = PurityContextFile.readWithQC(
                    mConfig.purpleQcFile(sampleId),
                    mConfig.purplePurityFile(sampleId)
            );

            dataItems.add(new DataItem(
                    DNA, ItemType.SAMPLE_TRAIT, GENDER.getAlias(),
                    boolToIntString(purityContext.gender() == Gender.MALE)
            ));

            dataItems.add(new DataItem(
                    DNA, ItemType.TUMOR_MUTATIONAL_BURDEN, MS_INDELS_TMB.getAlias(),
                    String.valueOf(purityContext.microsatelliteIndelsPerMb())
            ));

            dataItems.add(new DataItem(
                    DNA, ItemType.SAMPLE_TRAIT, WGD.getAlias(),
                    boolToIntString(purityContext.wholeGenomeDuplication())
            ));

            return dataItems;
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) sample traits - failed to load purity file from dir{}): {}",
                    sampleId,
                    mConfig.getPurpleDataDir(sampleId),
                    e.toString()
            );

            return null;
        }
    }
}
