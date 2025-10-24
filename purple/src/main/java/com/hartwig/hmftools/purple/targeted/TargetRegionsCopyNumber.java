package com.hartwig.hmftools.purple.targeted;

import java.util.List;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.TaggedRegion;

public record TargetRegionsCopyNumber(
        CobaltRatio cobaltRatio, List<TaggedRegion> overlappingRegions, PurpleCopyNumber purpleCopyNumber, GermlineStatus germlineStatus)
{


}
