package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isFrameshift;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isMissense;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isNonsense;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isSplice;

import java.util.Comparator;

import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.variant.SomaticVariant;

class TsgImpactComparator implements Comparator<SomaticVariant>
{
    @Override
    public int compare(final SomaticVariant o1, final SomaticVariant o2)
    {
        int firstWins = -1;
        int secondWins = 1;

        if(o1.type() == o2.type() && o1.canonicalCodingEffect() == o2.canonicalCodingEffect())
        {
            return 0;
        }

        if(DriverImpact.isFrameshift(o1))
        {
            return firstWins;
        }
        else if(DriverImpact.isFrameshift(o2))
        {
            return secondWins;
        }

        if(DriverImpact.isNonsense(o1))
        {
            return firstWins;
        }
        else if(DriverImpact.isNonsense(o2))
        {
            return secondWins;
        }

        if(DriverImpact.isSplice(o1))
        {
            return firstWins;
        }
        else if(DriverImpact.isSplice(o2))
        {
            return secondWins;
        }

        if(DriverImpact.isMissense(o1))
        {
            return firstWins;
        }
        else if(DriverImpact.isMissense(o2))
        {
            return secondWins;
        }

        throw new UnsupportedOperationException();
    }
}