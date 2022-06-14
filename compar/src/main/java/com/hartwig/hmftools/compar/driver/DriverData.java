package com.hartwig.hmftools.compar.driver;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class DriverData implements ComparableItem
{
    public final DriverCatalog DriverCatalog;

    protected static final String FLD_LIKELIHOOD = "Likelihood";
    protected static final String FLD_LIKE_METHOD = "LikelihoodMethod";

    public DriverData(final DriverCatalog driverCatalog)
    {
        DriverCatalog = driverCatalog;
    }

    @Override
    public Category category() { return DRIVER; }

    @Override
    public String key()
    {
        return format("%s_%s", DriverCatalog.driver(), DriverCatalog.gene());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", DriverCatalog.likelihoodMethod()));
        values.add(format("%.2f", DriverCatalog.driverLikelihood()));
        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DriverData otherDriver = (DriverData)other;

        if(!DriverCatalog.gene().equals(otherDriver.DriverCatalog.gene()))
            return false;

        if(DriverCatalog.driver() != otherDriver.DriverCatalog.driver())
        {
            if(DriverType.DRIVERS_LINX_SOMATIC.contains(DriverCatalog.driver())
            && DriverType.DRIVERS_LINX_SOMATIC.contains(otherDriver.DriverCatalog.driver()))
            {
                // matched due to 1.17 type expansion
            }
            else
            {
                return false;
            }
        }

        if(!DriverCatalog.transcript().equals(otherDriver.DriverCatalog.transcript()))
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final DriverData otherDriver = (DriverData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(
                diffs, FLD_LIKE_METHOD,
                DriverCatalog.likelihoodMethod().toString(), otherDriver.DriverCatalog.likelihoodMethod().toString());

        checkDiff(diffs, FLD_LIKELIHOOD, DriverCatalog.driverLikelihood(), otherDriver.DriverCatalog.driverLikelihood(), thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
