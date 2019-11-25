package com.hartwig.hmftools.sage.variant;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.sage.context.AltContext;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class SageVariantTierFactory {

    private final InPanel inPanel;
    private final InHotspot inHotspot;

    public SageVariantTierFactory(@NotNull final InPanel inPanel, @NotNull final InHotspot inHotspot) {
        this.inPanel = inPanel;
        this.inHotspot = inHotspot;
    }

    @NotNull
    public SageVariantTier tier(@NotNull final AltContext context) {
        if (inHotspot.isOnHotspot(context)) {
            return SageVariantTier.HOTSPOT;
        }
        if (inPanel.inPanel(context)) {
            return SageVariantTier.PANEL;
        }
        return SageVariantTier.WIDE;
    }

}
