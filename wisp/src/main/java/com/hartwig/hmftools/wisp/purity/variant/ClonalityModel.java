package com.hartwig.hmftools.wisp.purity.variant;

import java.util.List;

import com.hartwig.hmftools.wisp.common.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

public abstract class ClonalityModel
{
    protected final PurityConfig mConfig;
    protected final ResultsWriter mResultsWriter;

    protected final SampleData mSample;
    protected final List<SomaticVariant> mVariants;

    public ClonalityModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
        mVariants = variants;
    }

    abstract ClonalityResult calculate(final String sampleId, final FragmentCalcResult estimatedResult);

    public boolean useVariant(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return !variant.isFiltered() && !sampleFragData.isLowQual();
    }
}
