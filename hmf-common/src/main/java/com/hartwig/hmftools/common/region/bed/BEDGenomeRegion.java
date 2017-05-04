package com.hartwig.hmftools.common.region.bed;

import com.hartwig.hmftools.common.slicing.GenomeRegion;
import org.immutables.value.Value;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class BEDGenomeRegion implements GenomeRegion {
}
