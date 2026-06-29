package com.hartwig.hmftools.tars.liftback.tailextend;

import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MAX_EXTENSION;
import static com.hartwig.hmftools.tars.common.TarsConstants.DEFAULT_MIN_EXTENSION;

// Config for the tail-extension pass in TerminalReconciler. MinExtension is the minimum softclip length considered and
// the minimum bases reclaimed (one floor, since a reclaim can never exceed the softclip length). MaxExtension caps the
// walk to avoid consuming real unannotated junctions.
public class TailExtensionConfig
{
    public final boolean Enabled;
    public final int MinExtension;
    public final int MaxExtension;

    public TailExtensionConfig(final boolean enabled, final int minExtension, final int maxExtension)
    {
        Enabled = enabled;
        MinExtension = minExtension;
        MaxExtension = maxExtension;
    }

    public static TailExtensionConfig defaults()
    {
        return new TailExtensionConfig(false, DEFAULT_MIN_EXTENSION, DEFAULT_MAX_EXTENSION);
    }

    public static TailExtensionConfig enabledDefaults()
    {
        return new TailExtensionConfig(true, DEFAULT_MIN_EXTENSION, DEFAULT_MAX_EXTENSION);
    }
}
