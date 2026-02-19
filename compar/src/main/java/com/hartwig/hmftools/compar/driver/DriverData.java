package com.hartwig.hmftools.compar.driver;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class DriverData implements ComparableItem
{
    public final DriverCatalog DriverCatalog;
    public final String mComparisonChromosome;
    private final String mKey;
    private final boolean mCheckTranscript;

    protected static final String FLD_LIKELIHOOD = "Likelihood";
    protected static final String FLD_LIKE_METHOD = "LikelihoodMethod";
    protected static final String FLD_MIN_COPY_NUMBER = "MinCopyNumber";
    protected static final String FLD_MAX_COPY_NUMBER = "MaxCopyNumber";

    public DriverData(final DriverCatalog driverCatalog, final String comparisonChromosome, boolean checkTranscript)
    {
        DriverCatalog = driverCatalog;
        mComparisonChromosome = comparisonChromosome;
        mCheckTranscript = checkTranscript;

        String key = format("%s_%s", driverCatalog.driver(), driverCatalog.gene());
        mKey = driverCatalog.isCanonical() ? key : key + "_" + driverCatalog.transcript();
    }

    @Override
    public CategoryType category() { return DRIVER; }

    @Override
    public String key()
    {
        return mKey;
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", DriverCatalog.likelihoodMethod()));
        values.add(format("%.2f", DriverCatalog.driverLikelihood()));
        values.add(format("%.2f", DriverCatalog.minCopyNumber()));
        values.add(format("%.2f", DriverCatalog.maxCopyNumber()));
        return values;
    }

    @Override
    public boolean reportable() { return DriverCatalog.reportedStatus() == ReportedStatus.REPORTED; }

    @Override
    public String geneName() { return DriverCatalog.gene(); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DriverData otherDriver = (DriverData)other;

        if(!DriverCatalog.gene().equals(otherDriver.DriverCatalog.gene()))
            return false;

        if(DriverCatalog.driver() != otherDriver.DriverCatalog.driver())
            return false;

        if(mCheckTranscript || otherDriver.mCheckTranscript)
        {
            if(!DriverCatalog.transcript().equals(otherDriver.DriverCatalog.transcript())
            && DriverCatalog.isCanonical() != otherDriver.DriverCatalog.isCanonical())
            {
                return false;
            }
        }

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final DriverData otherDriver = (DriverData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(
                diffs, FLD_LIKE_METHOD,
                DriverCatalog.likelihoodMethod().toString(), otherDriver.DriverCatalog.likelihoodMethod().toString());

        checkDiff(diffs, FLD_LIKELIHOOD, DriverCatalog.driverLikelihood(), otherDriver.DriverCatalog.driverLikelihood(), thresholds);
        checkDiff(diffs, FLD_MIN_COPY_NUMBER, DriverCatalog.minCopyNumber(), otherDriver.DriverCatalog.minCopyNumber(), thresholds);
        checkDiff(diffs, FLD_MAX_COPY_NUMBER, DriverCatalog.maxCopyNumber(), otherDriver.DriverCatalog.maxCopyNumber(), thresholds);
        checkDiff(diffs, FLD_CHROMOSOME, mComparisonChromosome, otherDriver.mComparisonChromosome);
        checkDiff(diffs, FLD_CHROMOSOME_BAND, DriverCatalog.chromosomeBand(), otherDriver.DriverCatalog.chromosomeBand());

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
