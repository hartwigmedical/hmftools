package com.hartwig.hmftools.sage.context;

import java.util.Collections;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.primitives.Longs;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextMatch;

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

        int fullMatch = 0;
        int coreMatch = 0;
        for (ReadContextCandidate candidate : interimReadContexts) {
            final ReadContextMatch match = candidate.readContext().matchAtPosition(newReadContext);
            switch (match) {
                case FULL:
                    candidate.incrementFull(1);
                    fullMatch++;
                    break;
                case PARTIAL:
                case CORE:
                    candidate.incrementCore(1);
                    coreMatch++;
                    break;
            }
        }

        if (fullMatch == 0) {
            final ReadContextCandidate candidate = new ReadContextCandidate(newReadContext);
            candidate.incrementCore(coreMatch);
            interimReadContexts.add(candidate);
        }

    }

    public int readContextSupport() {
        return candidate.count();
    }

    public boolean finaliseAndValidate() {
        Collections.sort(interimReadContexts);
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

    static class ReadContextCandidate implements Comparable<ReadContextCandidate> {

        private final ReadContext readContext;
        private int fullMatch;
        private int coreMatch;

        ReadContextCandidate(@NotNull final ReadContext readContext) {
            assert (readContext.isComplete());
            this.readContext = readContext.minimiseFootprint();
        }

        public void incrementFull(int count) {
            fullMatch += count;
        }

        public void incrementCore(int count) {
            coreMatch += count;
        }

        public int count() {
            return fullMatch;
        }

        @NotNull
        public ReadContext readContext() {
            return readContext;
        }

        @Override
        public int compareTo(@NotNull final ReadContextCandidate o) {
            int fullCompare = -Integer.compare(fullMatch, o.fullMatch);
            if (fullCompare == 0) {
                return -Integer.compare(coreMatch, o.coreMatch);
            }

            return fullCompare;
        }
    }

}
