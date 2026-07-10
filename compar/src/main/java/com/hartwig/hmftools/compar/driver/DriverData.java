package com.hartwig.hmftools.compar.driver;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class DriverData implements ComparableItem
{
    public final DriverCatalog DriverCatalog;

    public final PurplePurity mPurity;
    public final String mComparisonChromosome;
    private final String mKey;
    private final boolean mCheckTranscript;

    public DriverData(
            final DriverCatalog driverCatalog, final PurplePurity purity, final String comparisonChromosome, boolean checkTranscript)
    {
        DriverCatalog = driverCatalog;

        mPurity = purity;
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
}
