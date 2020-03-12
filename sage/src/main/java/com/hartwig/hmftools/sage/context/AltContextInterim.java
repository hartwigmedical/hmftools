package com.hartwig.hmftools.sage.context;

import java.util.Comparator;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.primitives.Longs;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

public class AltContextInterim implements VariantHotspot {

    private final RefContext refContext;
    private final String ref;
    private final String alt;
    private final List<ReadContextCounter> interimReadContexts = Lists.newArrayList();

    private ReadContextCounter readContextCounter;
    private int rawSupportAlt;
    private int rawBaseQualityAlt;

    public AltContextInterim(final RefContext refContext, final String ref, final String alt) {
        this.refContext = refContext;
        this.ref = ref;
        this.alt = alt;
    }

    public void incrementAltRead(int baseQuality) {
        this.rawSupportAlt++;
        this.rawBaseQualityAlt += baseQuality;
    }

    public void addReadContext(@NotNull final ReadContext newReadContext) {
        if (readContextCounter != null) {
            throw new IllegalStateException();
        }

        if (newReadContext.isComplete()) {
            boolean readContextMatch = false;
            for (ReadContextCounter counter : interimReadContexts) {
                if (counter.incrementCounters(newReadContext)) {
                    readContextMatch = true;
                    break;
                }

            }

            if (!readContextMatch) {
                interimReadContexts.add(new ReadContextCounter(this, newReadContext));
            }
        }
    }

    @Nullable
    public ReadContextCounter primaryReadContext() {
        interimReadContexts.sort(Comparator.comparingInt(ReadContextCounter::altSupport).reversed());
        return interimReadContexts.isEmpty()
                ? null
                : new ReadContextCounter(this, interimReadContexts.get(0).readContext().minimiseFootprint());
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

    public int rawSupportRef() {
        return refContext.rawSupportRef();
    }

    public int rawBaseQualityRef() {
        return refContext.rawBaseQualityRef();
    }

    public int rawSupportAlt() {
        return rawSupportAlt;
    }

    public int rawDepth() {
        return refContext.rawDepth();
    }

    public double rawVaf() {
        return refContext.rawDepth() == 0 ? 0 : ((double) rawSupportAlt) / rawDepth();
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
