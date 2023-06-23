package com.hartwig.hmftools.ctdna.purity.variant;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;
import com.hartwig.hmftools.ctdna.purity.ResultsWriter;
import com.hartwig.hmftools.ctdna.purity.SampleData;

public class DropoutRateModel
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final SampleData mSample;
    private final List<SomaticVariant> mVariants;

    public DropoutRateModel(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
        mVariants = Lists.newArrayList();
    }



}
