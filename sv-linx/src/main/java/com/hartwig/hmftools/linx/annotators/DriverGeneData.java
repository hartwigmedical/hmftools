package com.hartwig.hmftools.linx.annotators;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.AMP;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DEL;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

public class DriverGeneData
{
    public final DriverCatalog DriverData;
    public final EnsemblGeneData GeneData;
    public final TranscriptData TransData;
    public final GeneCopyNumber GeneCN;

    private final List<DriverGeneEvent> mEvents;

    public DriverGeneData(final DriverCatalog driverData, final EnsemblGeneData geneData,
            final TranscriptData transData, final GeneCopyNumber geneCN)
    {
        DriverData = driverData;
        GeneData = geneData;
        TransData = transData;
        GeneCN = geneCN;

        mEvents = Lists.newArrayList();
    }

    public final List<DriverGeneEvent> getEvents() {return mEvents; }
    public void addEvent(final DriverGeneEvent event) { mEvents.add(event); }

    public boolean fullyMatched()
    {
        if(DriverData.driver() == AMP)
            return !mEvents.isEmpty();

        if(DriverData.driver() == DEL)
        {
            if(DriverData.likelihoodMethod() == LikelihoodMethod.DEL)
            {
                return mEvents.size() >= 2;
            }
            else
            {
                return !mEvents.isEmpty();
            }
        }

        return false;
    }

    public String toString()
    {
        return String.format("%s: %s ",
            DriverData.likelihoodMethod(), DriverData.gene());
    }
}
