package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

import com.hartwig.hmftools.sage.common.SimpleVariant;

import org.jetbrains.annotations.NotNull;

public class HotspotSelector
{
    private final PositionSelector<SimpleVariant> mHotspotSelector;

    public HotspotSelector(final List<SimpleVariant> hotspots)
    {
        mHotspotSelector = new PositionSelector<>(hotspots);
    }

    public boolean isHotspot(final SimpleVariant variant)
    {
        final AtomicBoolean hotspotMatch = new AtomicBoolean(false);

        mHotspotSelector.select(variant.position(), variant.position(), hotspot ->
        {
            if(hotspot.alt().equals(variant.alt()) && hotspot.ref().equals(variant.ref()))
            {
                hotspotMatch.set(true);
            }
        });

        return hotspotMatch.get();
    }
}
