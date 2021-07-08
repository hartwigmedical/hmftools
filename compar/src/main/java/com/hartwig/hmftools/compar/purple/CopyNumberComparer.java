package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.DRIVER;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.compar.driver.DriverData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.compress.utils.Lists;

public class CopyNumberComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public CopyNumberComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = mConfig.Categories.get(DRIVER);

        final List<List<ComparableItem>> sourceDrivers = Lists.newArrayList();

        for(String sourceName : mConfig.DbSourceNames)
        {
            sourceDrivers.add(getSampleDrivers(sampleId, mConfig.DbConnections.get(sourceName)));
        }

        for(int i = 0; i < mConfig.DbSourceNames.size() - 1; ++i)
        {
            final String source1 = mConfig.DbSourceNames.get(i);

            for(int j = i + 1; j < mConfig.DbSourceNames.size(); ++j)
            {
                final String source2 = mConfig.DbSourceNames.get(j);

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, sourceDrivers.get(i), sourceDrivers.get(j));
            }
        }
    }

    private List<ComparableItem> getSampleDrivers(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<DriverCatalog> drivers = dbAccess.readDriverCatalog(sampleId);
        final List<LinxDriver> svDrivers = dbAccess.readSvDriver(sampleId);

        final List<ComparableItem> driverDataList = Lists.newArrayList();

        for(DriverCatalog driver : drivers)
        {
            List<LinxDriver> svDriverList = svDrivers.stream().filter(x -> x.gene().equals(driver.gene())).collect(Collectors.toList());
            driverDataList.add(new DriverData(driver, svDriverList));
        }

        return driverDataList;
    }


}
