package com.hartwig.hmftools.compar.resource;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.compar.common.MismatchType;

public record ResourceMismatch(MismatchType type, String key, List<String> diffValues, List<String> oldLines, List<String> newLines)
{
    public ResourceMismatch(final MismatchType type)
    {
        this(type, "", Collections.emptyList(), Collections.emptyList(), Collections.emptyList());
    }
}
