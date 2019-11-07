package com.hartwig.hmftools.sage.context;

import java.util.Comparator;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.primitives.Longs;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
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
    private int altReads;


    public AltContext(final String sample, final VariantHotspot hotspot) {
        refContext = new RefContext(sample, hotspot.chromosome(), hotspot.position());
        this.alt = hotspot.alt();
        this.ref = hotspot.ref();
    }

    public AltContext(final RefContext refContext, final String ref, final String alt) {
        this.refContext = refContext;
        this.ref = ref;
        this.alt = alt;
    }

    public RefContext refContext() {
        return refContext;
    }

    public void incrementAltRead() {
        this.altReads++;
    }

    public void addReadContext(@NotNull final ReadContext newReadContext) {
        if (readContextCounter != null) {
            throw new IllegalStateException();
        }

        if (newReadContext.isComplete()) {
            boolean readContextMatch = false;
            for (ReadContextCounter counter : interimReadContexts) {
                readContextMatch |= counter.incrementCounters(newReadContext);
            }

            if (!readContextMatch) {
                interimReadContexts.add(new ReadContextCounter(this, newReadContext));
            }
        }
    }

    @NotNull
    public ReadContextCounter setPrimaryReadCounterFromInterim() {
        interimReadContexts.sort(Comparator.comparingInt(ReadContextCounter::support).reversed());
        readContextCounter = interimReadContexts.isEmpty()
                ? new ReadContextCounter(this, ReadContextFactory.dummy((int) position(), alt))
                : new ReadContextCounter(this, interimReadContexts.get(0).readContext());

        interimReadContexts.clear();

        return primaryReadContext();
    }

    @NotNull
    public ReadContextCounter primaryReadContext() {
        if (readContextCounter == null) {
            throw new IllegalStateException();
        }

        return readContextCounter;
    }

    public void setPrimaryReadContext(@NotNull final ReadContextCounter readContext) {
        readContextCounter = readContext;
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

    public int altSupport() {
        return altReads;
    }

    public int refSupport() {
        return refContext.refReads();
    }

    public int readDepth() {
        return refContext.readDepth();
    }

    public int subprimeReadDepth() {
        return refContext.subprimeReadDepth();
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
