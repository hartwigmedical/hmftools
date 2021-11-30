package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.gripss.filters.FilterConstants.SHORT_RESCUE_LENGTH;

public class LinkRescue
{
    public static boolean tooShortToRescue(int length) { return length < SHORT_RESCUE_LENGTH; }
}

