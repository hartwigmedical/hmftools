package com.hartwig.hmftools.esvee;

import org.immutables.value.Value;

@Value.Immutable
public abstract class Junction {
    public abstract String chromosome();
    public abstract int position();
    public abstract Direction orientation();
    public abstract int junctionFragments();
    public abstract int supportFragments();
    public abstract int discordantFragments();
    public abstract int lowMapQualityFragments();
    public abstract int maxMapQuality();
    public abstract int maxSoftClipLength();
    public abstract int baseDepth();
    public abstract boolean hasPolyAT();
    public abstract boolean isIndel();
    public abstract boolean isHotspot();
    public abstract String softClipBases();
    public abstract String initialReadId();

    @Override
    public String toString()
    {
        return String.format("%s:%s (%s)", chromosome(), position(), orientation());
    }
}
