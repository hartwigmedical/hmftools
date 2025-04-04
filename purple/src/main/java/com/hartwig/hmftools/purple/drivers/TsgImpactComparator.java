package com.hartwig.hmftools.purple.drivers;

import java.util.Comparator;

import com.hartwig.hmftools.common.driver.DriverImpact;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public class TsgImpactComparator implements Comparator<SomaticVariant>
{
    @Override
    public int compare(final SomaticVariant o1, final SomaticVariant o2)
    {
        int firstWins = -1;
        int secondWins = 1;

        CodingEffect codingEffect1 = o1.variantImpact().CanonicalCodingEffect;
        CodingEffect codingEffect2 = o2.variantImpact().CanonicalCodingEffect;

        DriverImpact impact1 = DriverImpact.select(o1.type(), codingEffect1);
        DriverImpact impact2 = DriverImpact.select(o2.type(), codingEffect2);

        int impactScore1 = driverImpactScore(impact1);
        int impactScore2 = driverImpactScore(impact2);

        if(impactScore1 < impactScore2)
            return firstWins;
        else if(impactScore1 > impactScore2)
            return secondWins;

        if(o1.position() == o2.position())
            return 0;

        return o1.position() < o2.position() ? firstWins : secondWins;
    }

    public static int driverImpactScore(final DriverImpact impact)
    {
        switch(impact)
        {
            case FRAMESHIFT: return 0;
            case NONSENSE: return 1;
            case SPLICE: return 2;
            case MISSENSE: return 3;
            case INFRAME: return 4;
        }

        return 5;
    }
}