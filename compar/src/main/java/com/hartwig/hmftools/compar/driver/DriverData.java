package com.hartwig.hmftools.compar.driver;

import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.CommonUtils.diffValue;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;

import org.apache.commons.compress.utils.Lists;

public class DriverData implements ComparableItem
{
    public final DriverCatalog DriverCatalog;

    public DriverData(final DriverCatalog driverCatalog)
    {
        DriverCatalog = driverCatalog;
    }

    public Category category() { return DRIVER; }

    public boolean reportable() { return true; }

    public boolean matches(final ComparableItem other)
    {
        final DriverData otherDriver = (DriverData)other;

        if(DriverCatalog.driver() != otherDriver.DriverCatalog.driver())
            return false;

        if(!DriverCatalog.gene().equals(otherDriver.DriverCatalog.gene()))
            return false;

        return true;

    }

    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final DriverData otherDriver = (DriverData)other;

        final List<String> diffs = Lists.newArrayList();

        if(matchLevel == REPORTABLE)
            return diffs;

        if(diffValue(DriverCatalog.driverLikelihood(), otherDriver.DriverCatalog.driverLikelihood()))
        {
            diffs.add(String.format("likelihood(%.3f/%.3f)",
                    DriverCatalog.driverLikelihood(), otherDriver.DriverCatalog.driverLikelihood()));
        }

        return diffs;
    }

    public String description()
    {
        return String.format("%s_%s", DriverCatalog.driver(), DriverCatalog.gene());
    }
}
