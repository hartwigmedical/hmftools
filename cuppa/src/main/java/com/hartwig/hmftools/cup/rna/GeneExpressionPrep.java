package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.PrepConfig;
import com.hartwig.hmftools.cup.traits.SampleTraitsData;

public class GeneExpressionPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    public GeneExpressionPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return CategoryType.GENE_EXP; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        final String purpleDataDir = formSamplePath(mConfig.PurpleDir, sampleId);

        try
        {
            final PurityContext purityContext = PurityContextFile.read(purpleDataDir, sampleId);

            SampleTraitsData traitsData = SampleTraitsData.from(sampleId, purityContext, 0);


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
