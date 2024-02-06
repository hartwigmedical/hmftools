package com.hartwig.hmftools.cup.traits;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupCalcs.boolToIntString;
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
import com.hartwig.hmftools.cup.common.CupCalcs;
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

    private static final String GENDER_IS_MALE = "IS_MALE";

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        final String purpleDataDir = mConfig.getPurpleDataDir(sampleId);

        try
        {
            final PurityContext purityContext = PurityContextFile.read(purpleDataDir, sampleId);

            dataItems.add(new DataItem(DNA, ItemType.SAMPLE_TRAIT, GENDER_IS_MALE, boolToIntString(purityContext.gender() == Gender.MALE)));

            dataItems.add(new DataItem(DNA, ItemType.SAMPLE_TRAIT, MS_INDELS_TMB.toString(), String.valueOf(purityContext.microsatelliteIndelsPerMb())));

            dataItems.add(new DataItem(DNA, ItemType.SAMPLE_TRAIT, WGD.toString(), boolToIntString(purityContext.wholeGenomeDuplication())));

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
