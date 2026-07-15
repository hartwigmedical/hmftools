package com.hartwig.hmftools.compar.driver;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.CommonUtils.findDiffs;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class DriverData implements ComparableItem
{
    protected static final String FLD_LIKELIHOOD = "Likelihood";
    protected static final String FLD_LIKE_METHOD = "LikelihoodMethod";
    protected static final String FLD_MIN_COPY_NUMBER = "MinCopyNumber";
    protected static final String FLD_MAX_COPY_NUMBER = "MaxCopyNumber";

    public final DriverCatalog DriverCatalog;

    public final PurplePurity mPurity;
    public final String mComparisonChromosome;
    private final String mKey;
    private final boolean mCheckTranscript;
    private final boolean mIsPass;

    public DriverData(final DriverCatalog driverCatalog, final PurplePurity purity, final String comparisonChromosome,
            boolean checkTranscript, boolean isPass)
    {
        DriverCatalog = driverCatalog;

        mPurity = purity;
        mComparisonChromosome = comparisonChromosome;
        mCheckTranscript = checkTranscript;

        String key = format("%s_%s", driverCatalog.driver(), driverCatalog.gene());
        mKey = driverCatalog.isCanonical() ? key : key + "_" + driverCatalog.transcript();
        mIsPass = isPass;
    }

    @Override
    public CategoryType category() { return DRIVER; }

    @Override
    public String key()
    {
        return mKey;
    }

    @Override
    public boolean isPass()
    {
        return mIsPass;
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

    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final FieldConfig fieldConfig,
            final boolean includeMatches)
    {
        final DriverData otherData = (DriverData) other;
        final List<String> diffs = Lists.newArrayList();
        List<String> alwaysCompareFields = List.of(
                FLD_LIKE_METHOD, FLD_MIN_COPY_NUMBER, FLD_MAX_COPY_NUMBER, FLD_CHROMOSOME, FLD_CHROMOSOME_BAND);
        diffs.addAll(findDiffs(this, otherData, fieldConfig.getFields(category(), alwaysCompareFields)));

        if(isPass() && otherData.isPass())
        {
            diffs.addAll(findDiffs(this, otherData, fieldConfig.getFields(category(), List.of(FLD_LIKELIHOOD))));
        }
        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
