package com.hartwig.hmftools.purple.drivers;

import java.util.Comparator;

import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.purple.somatic.SomaticData;

public class TsgImpactComparator implements Comparator<SomaticData>
{
    @Override
    public int compare(final SomaticData o1, final SomaticData o2)
    {
        int firstWins = -1;
        int secondWins = 1;

        CodingEffect codingEffect1 = o1.variantImpact().CanonicalCodingEffect;
        CodingEffect codingEffect2 = o2.variantImpact().CanonicalCodingEffect;

        if(codingEffect1 == codingEffect2 && o1.type() == o2.type())
            return 0;

        if(DriverImpact.isFrameshift(o1.type(), codingEffect1))
        {
            return firstWins;
        }
        else if(DriverImpact.isFrameshift(o2.type(), codingEffect2))
        {
            return secondWins;
        }

        if(DriverImpact.isNonsense(o1.type(), codingEffect1))
        {
            return firstWins;
        }
        else if(DriverImpact.isNonsense(o2.type(), codingEffect2))
        {
            return secondWins;
        }

        if(DriverImpact.isSplice(codingEffect1))
        {
            return firstWins;
        }
        else if(DriverImpact.isSplice(codingEffect2))
        {
            return secondWins;
        }

        if(DriverImpact.isMissense(o1.type(), codingEffect1))
        {
            return firstWins;
        }
        else if(DriverImpact.isMissense(o2.type(), codingEffect2))
        {
            return secondWins;
        }

        if(codingEffect1 == codingEffect2)
            return 0;

        throw new UnsupportedOperationException();
    }
}