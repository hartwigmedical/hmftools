package com.hartwig.hmftools.wisp.variant;

import java.util.List;

import com.hartwig.hmftools.wisp.SampleData;
import com.hartwig.hmftools.wisp.WispConfig;
import com.hartwig.hmftools.wisp.ResultsWriter;

public abstract class ClonalityModel
{
    protected final WispConfig mConfig;
    protected final ResultsWriter mResultsWriter;

    protected final SampleData mSample;
    protected final List<SomaticVariant> mVariants;

    public ClonalityModel(
            final WispConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
        mVariants = variants;
    }

    abstract void calculate(final String sampleId, final FragmentTotals fragmentTotals, final PurityCalcData purityCalcData, double noiseRate);

    public boolean useVariant(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return !variant.isFiltered() && !sampleFragData.isFiltered();
    }
}
