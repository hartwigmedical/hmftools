package com.hartwig.hmftools.ctdna.purity.variant;

import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MIN_QUAL_PER_AD;

import java.util.List;

import com.hartwig.hmftools.ctdna.common.SampleData;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;
import com.hartwig.hmftools.ctdna.purity.ResultsWriter;

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
        return variant.PassFilters
                && variant.sequenceGcRatio() >= mConfig.GcRatioMin
                && (sampleFragData.qualPerAlleleFragment() > MIN_QUAL_PER_AD || sampleFragData.AlleleCount == 0);
    }
}
