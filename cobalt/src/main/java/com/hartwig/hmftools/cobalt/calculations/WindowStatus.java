package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public final class WindowStatus extends ChrBaseRegion
{
    private final boolean Excluded;
    private final boolean Unmappable;

    public WindowStatus(GCProfile gcProfile, boolean excluded)
    {
        super(gcProfile.chromosome(), gcProfile.start(), gcProfile.end());
        Excluded = excluded;
        Unmappable = !gcProfile.isMappable();
    }

    public boolean maskedOut()
    {
        return Excluded || Unmappable;
    }
}
