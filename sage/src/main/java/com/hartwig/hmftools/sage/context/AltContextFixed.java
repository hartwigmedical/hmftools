package com.hartwig.hmftools.sage.context;

import javax.annotation.Nullable;

import com.google.common.primitives.Longs;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

public class AltContextFixed implements AltContext {

    private final String ref;
    private final String alt;
    private final RefContext refContext;
    private final ReadContextCounter readContextCounter;
    private final boolean realign;

    private int rawSupportAlt;
    private int rawBaseQualityAlt;

    public AltContextFixed(final RefContext refContext, final String ref, final String alt, boolean realign,
            final ReadContext readContext) {
        this.refContext = refContext;
        this.ref = ref;
        this.alt = alt;
        this.readContextCounter = new ReadContextCounter(refContext.sample(), this, readContext, realign);
        this.realign = realign;
    }

    public boolean realign() {
        return realign;
    }

    public void incrementAltRead(int baseQuality) {
        this.rawSupportAlt++;
        this.rawBaseQualityAlt += baseQuality;
    }

    @NotNull
    public ReadContextCounter primaryReadContext() {
        return readContextCounter;
    }

    @NotNull
    @Override
    public String ref() {
        return ref;
    }

    @NotNull
    @Override
    public String alt() {
        return alt;
    }

    @NotNull
    @Override
    public String chromosome() {
        return refContext.chromosome();
    }

    @Override
    public long position() {
        return refContext.position();
    }

    public int rawSupportAlt() {
        return rawSupportAlt;
    }

    public int rawDepth() {
        return refContext.rawDepth();
    }

    public int rawBaseQualityAlt() {
        return rawBaseQualityAlt;
    }

    @NotNull
    public String sample() {
        return refContext.sample();
    }

    @Override
    public boolean equals(@Nullable Object another) {
        if (this == another) {
            return true;
        }
        return another instanceof VariantHotspot && equalTo((VariantHotspot) another);
    }

    private boolean equalTo(VariantHotspot another) {
        return ref().equals(another.ref()) && alt().equals(another.alt()) && chromosome().equals(another.chromosome())
                && position() == another.position();
    }

    @Override
    public int hashCode() {
        int h = 5381;
        h += (h << 5) + ref.hashCode();
        h += (h << 5) + alt.hashCode();
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Longs.hashCode(position());
        return h;
    }
}
