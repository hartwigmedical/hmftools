package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public final class WindowStatus extends ChrBaseRegion
{
    private final boolean Excluded;
    private final boolean Unmappable;
    private final boolean NonDiploid;

    public WindowStatus(GCProfile gcProfile, boolean excluded, boolean nonDiploid)
    {
        super(gcProfile.chromosome(), gcProfile.start(), gcProfile.end());
        Excluded = excluded;
        Unmappable = !gcProfile.isMappable();
        this.NonDiploid = nonDiploid;
    }

    public boolean maskedOut()
    {
        return Excluded || Unmappable || NonDiploid;
    }

    @Override
    public String toString()
    {
        return String.format("%s:%d-%d %b %b %b", Chromosome, start(), end(), Excluded, Unmappable, NonDiploid);

    }
}
