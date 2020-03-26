package com.hartwig.hmftools.sage.context;

import java.util.Comparator;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.primitives.Longs;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextFactory;

import org.jetbrains.annotations.NotNull;

public class AltContext implements VariantHotspot {

    private final RefContext refContext;
    private final String ref;
    private final String alt;
    private final List<ReadContextCounter> interimReadContexts = Lists.newArrayList();

    private ReadContextCounter readContextCounter;
    private int rawSupportAlt;
    private int rawBaseQualityAlt;

    public AltContext(final RefContext refContext, final String ref, final String alt) {
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
                interimReadContexts.add(new ReadContextCounter(refContext.sample(), this, newReadContext, false));
            }
        }
    }

    @NotNull
    public ReadContextCounter setPrimaryReadCounterFromInterim() {
        interimReadContexts.sort(Comparator.comparingInt(ReadContextCounter::altSupport).reversed());
        readContextCounter = interimReadContexts.isEmpty()
                ? new ReadContextCounter(refContext.sample(), this, ReadContextFactory.dummy((int) position(), alt), false)
                : new ReadContextCounter(refContext.sample(), this, interimReadContexts.get(0).readContext().minimiseFootprint(), false);

        interimReadContexts.clear();

        return primaryReadContext();
    }

    @NotNull
    public ReadContextCounter primaryReadContext() {
        if (readContextCounter == null) {
            setPrimaryReadCounterFromInterim();
        }

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

    public int rawAltSupport() {
        return rawSupportAlt;
    }

    public int rawDepth() {
        return refContext.rawDepth();
    }


    public int rawAltBaseQuality() {
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
