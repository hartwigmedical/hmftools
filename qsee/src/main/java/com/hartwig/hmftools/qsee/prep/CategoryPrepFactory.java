package com.hartwig.hmftools.qsee.prep;

import java.util.List;

import com.hartwig.hmftools.qsee.prep.category.GcBiasPrep;
import com.hartwig.hmftools.qsee.prep.category.CoverageDistributionPrep;
import com.hartwig.hmftools.qsee.prep.category.FragLengthDistributionPrep;
import com.hartwig.hmftools.qsee.prep.category.MissedGeneVariantPrep;
import com.hartwig.hmftools.qsee.prep.category.BaseQualRecalibrationPrep;
import com.hartwig.hmftools.qsee.prep.category.DuplicateFreqPrep;
import com.hartwig.hmftools.qsee.prep.category.MsIndelErrorPrep;
import com.hartwig.hmftools.qsee.prep.category.SummaryTablePrep;

public class CategoryPrepFactory
{
    CommonPrepConfig mConfig;

    public CategoryPrepFactory(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public List<CategoryPrep> createCategoryPreps()
    {
        return List.of(
                new SummaryTablePrep(mConfig),
                new CoverageDistributionPrep(mConfig),
                new FragLengthDistributionPrep(mConfig),
                new MissedGeneVariantPrep(mConfig),
                new DuplicateFreqPrep(mConfig),
                new GcBiasPrep(mConfig),
                new BaseQualRecalibrationPrep(mConfig),
                new MsIndelErrorPrep(mConfig)
        );
    }
}
