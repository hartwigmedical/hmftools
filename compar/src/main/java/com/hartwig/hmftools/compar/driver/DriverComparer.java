package com.hartwig.hmftools.compar.driver;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DriverComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public DriverComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        for(Map.Entry<String, DatabaseAccess> entry1 : mConfig.DbConnections.entrySet())
        {
            final String source1 = entry1.getKey();
            final List<ComparableItem> drivers1 = getSampleDrivers(sampleId, entry1.getValue());

            for(Map.Entry<String, DatabaseAccess> entry2 : mConfig.DbConnections.entrySet())
            {
                if(entry1 == entry2)
                    continue;

                final String source2 = entry2.getKey();

                final List<ComparableItem> drivers2 = getSampleDrivers(sampleId, entry2.getValue());

                CommonUtils.compareItems(sampleId, mismatches, source1, source2, drivers1, drivers2);
            }
        }
    }

    private List<ComparableItem> getSampleDrivers(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<DriverCatalog> drivers = dbAccess.readDriverCatalog(sampleId);

        return drivers.stream().map(x -> new DriverData(x)).collect(Collectors.toList());
    }


}
