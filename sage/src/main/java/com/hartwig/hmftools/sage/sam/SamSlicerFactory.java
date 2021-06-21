package com.hartwig.hmftools.sage.sam;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

public class SamSlicerFactory
{
    private final SageConfig mConfig;
    private final List<BaseRegion> mPanel;

    public SamSlicerFactory(@NotNull final SageConfig config, @NotNull final List<BaseRegion> panel)
    {
        mConfig = config;
        mPanel = panel;
    }

    @NotNull
    public SamSlicer create(final BaseRegion slice)
    {
        return mConfig.PanelOnly ? panelOnly(slice) : fullSlice(slice);
    }

    @NotNull
    private SamSlicer fullSlice(@NotNull final BaseRegion slice)
    {
        return new SamSlicer(0, slice);
    }

    @NotNull
    private SamSlicer panelOnly(@NotNull final BaseRegion slice)
    {
        return new SamSlicer(0, slice, mPanel);
    }
}
