package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.MismatchType.PRESENCE;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DriverComparator implements Comparator
{
    private final ComparConfig mConfig;

    DriverComparator(final ComparConfig config)
    {
        mConfig = config;

    }

    public void processSample(final String sampleId, final List<DataMismatch> mismatches)
    {
        for(Map.Entry<String, DatabaseAccess> entry1 : mConfig.DbConnections.entrySet())
        {
            for(Map.Entry<String, DatabaseAccess> entry2 : mConfig.DbConnections.entrySet())
            {
                if(entry1 == entry2)
                    continue;

                final String source1 = entry1.getKey();
                final String source2 = entry2.getKey();

                final List<DriverCatalog> drivers1 = entry1.getValue().readDriverCatalog(sampleId);
                final List<DriverCatalog> drivers2 = entry2.getValue().readDriverCatalog(sampleId);

                int index1 = 0;
                while(index1 < drivers1.size())
                {
                    final DriverCatalog driver1 = drivers1.get(index1);

                    boolean matched = false;

                    int index2 = 0;
                    while(index2 < drivers2.size())
                    {
                        final DriverCatalog driver2 = drivers2.get(index2);

                        if(matched(driver1, driver2))
                        {
                            drivers1.remove(index1);
                            drivers2.remove(index2);
                            matched = true;
                            break;
                        }
                        else
                        {
                            ++index2;
                        }
                    }

                    if(!matched)
                        ++index1;
                }

                for(DriverCatalog driver : drivers1)
                {
                    mismatches.add(new DataMismatch(
                            sampleId, DRIVER, PRESENCE, source1, source2, String.format("%s_%s", driver.driver(), driver.gene()), ""));
                }

                for(DriverCatalog driver : drivers2)
                {
                    mismatches.add(new DataMismatch(
                            sampleId, DRIVER, PRESENCE, source2, source1, String.format("%s_%s", driver.driver(), driver.gene()), ""));
                }
            }
        }
    }

    private boolean matched(final DriverCatalog driver1, final DriverCatalog driver2)
    {
        if(driver1.driver() != driver2.driver())
            return false;

        if(!driver1.gene().equals(driver2.gene()))
            return false;

        return true;
    }


}
