package com.hartwig.hmftools.sage.context;

import java.util.Comparator;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.primitives.Longs;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class AltContext implements VariantHotspot {

    private final String ref;
    private final String alt;
    private final RefContext refContext;
    private final List<ReadContextCandidate> interimReadContexts = Lists.newArrayList();

    private int rawSupportAlt;
    private int rawBaseQualityAlt;
    private ReadContextCandidate candidate;

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
        if (candidate != null) {
            throw new IllegalStateException();
        }

        if (newReadContext.isComplete()) {
            boolean readContextMatch = false;
            for (ReadContextCandidate counter : interimReadContexts) {
                if (counter.incrementCounters(newReadContext)) {
                    readContextMatch = true;
                    break;
                }
            }

            if (!readContextMatch) {
                interimReadContexts.add(new ReadContextCandidate(newReadContext));
            }
        }
    }

    public int readContextSupport() {
        return candidate.count();
    }

    public boolean finaliseAndValidate() {
        interimReadContexts.sort(Comparator.comparingInt(ReadContextCandidate::count).reversed());
        if (!interimReadContexts.isEmpty()) {
            candidate = interimReadContexts.get(0);
        }
        interimReadContexts.clear();
        return candidate != null && candidate.readContext().isComplete();
    }

    public ReadContext readContext() {
        return candidate.readContext();
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

    static class ReadContextCandidate {

        private final ReadContext readContext;
        private int count;

        ReadContextCandidate(@NotNull final ReadContext readContext) {
            assert (readContext.isComplete());
            this.readContext = readContext;
        }

        public boolean incrementCounters(@NotNull final ReadContext other) {
            if (readContext.isFullMatch(other)) {
                count++;
                return true;
            }

            return false;
        }

        public int count() {
            return count;
        }

        @NotNull
        public ReadContext readContext() {
            return readContext;
        }
    }

}
