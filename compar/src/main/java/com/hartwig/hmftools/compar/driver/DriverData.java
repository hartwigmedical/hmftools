package com.hartwig.hmftools.compar.driver;

import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class DriverData implements ComparableItem
{
    public final DriverCatalog DriverCatalog;
    public final String MappedGeneName; // overridden with new mapped name if applicable
    public final List<LinxDriver> SvDrivers;

    public DriverData(final DriverCatalog driverCatalog, final List<LinxDriver> svDrivers, final String mappedName)
    {
        DriverCatalog = driverCatalog;
        SvDrivers = svDrivers;
        MappedGeneName = mappedName;
    }

    @Override
    public Category category() { return DRIVER; }

    @Override
    public String key()
    {
        return String.format("%s_%s_%s", DriverCatalog.driver(), DriverCatalog.gene(), DriverCatalog.likelihoodMethod());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DriverData otherDriver = (DriverData)other;

        //if(!DriverCatalog.gene().equals(otherDriver.DriverCatalog.gene()))

        if(!MappedGeneName.equals(otherDriver.MappedGeneName))
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


        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel)
    {
        return null;
    }

    @Override
    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final DriverData otherDriver = (DriverData)other;

        final List<String> diffs = Lists.newArrayList();

        if(matchLevel == REPORTABLE)
            return diffs;

        checkDiff(
                diffs, "likelihoodMethod",
                DriverCatalog.likelihoodMethod().toString(), otherDriver.DriverCatalog.likelihoodMethod().toString());

        checkDiff(diffs, "likelihood", DriverCatalog.driverLikelihood(), otherDriver.DriverCatalog.driverLikelihood());

        // check matches in Linx cluster event types
        boolean hasDiffs = ((DriverData) other).SvDrivers.size() != SvDrivers.size();

        if(!hasDiffs)
        {
            for(final LinxDriver svDriver : SvDrivers)
            {
                if(otherDriver.SvDrivers.stream().noneMatch(x -> x.eventType().equals(svDriver.eventType())))
                {
                    hasDiffs = true;
                    break;
                }
            }
        }

        if(hasDiffs)
        {
            final StringJoiner eventTypes = new StringJoiner(ITEM_DELIM);
            SvDrivers.stream().map(x -> x.eventType()).forEach(x -> eventTypes.add(x));

            final StringJoiner otherEventTypes = new StringJoiner(ITEM_DELIM);
            otherDriver.SvDrivers.stream().map(x -> x.eventType()).forEach(x -> otherEventTypes.add(x));

            diffs.add(String.format("svDriver eventTypes(%s/%s)",
                    eventTypes.toString(), otherEventTypes.toString()));
        }


        return diffs;
    }

}
