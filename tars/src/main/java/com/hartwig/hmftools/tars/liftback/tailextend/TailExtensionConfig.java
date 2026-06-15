package com.hartwig.hmftools.tars.liftback.tailextend;

// Config for SoftclipTailExtender. MaxExtension caps the walk to avoid consuming real unannotated junctions.
public class TailExtensionConfig
{
    public static final int DEFAULT_MIN_SOFTCLIP_LENGTH = 3;
    public static final int DEFAULT_MIN_EXTENSION = 3;
    public static final int DEFAULT_MAX_EXTENSION = 30;

    public final boolean Enabled;
    public final int MinSoftclipLength;
    public final int MinExtension;
    public final int MaxExtension;

    public TailExtensionConfig(
            final boolean enabled, final int minSoftclipLength,
            final int minExtension, final int maxExtension)
    {
        Enabled = enabled;
        MinSoftclipLength = minSoftclipLength;
        MinExtension = minExtension;
        MaxExtension = maxExtension;
    }

    public static TailExtensionConfig defaults()
    {
        return new TailExtensionConfig(
                false, DEFAULT_MIN_SOFTCLIP_LENGTH, DEFAULT_MIN_EXTENSION, DEFAULT_MAX_EXTENSION);
    }

    public static TailExtensionConfig enabledDefaults()
    {
        return new TailExtensionConfig(
                true, DEFAULT_MIN_SOFTCLIP_LENGTH, DEFAULT_MIN_EXTENSION, DEFAULT_MAX_EXTENSION);
    }
}
