package com.hartwig.hmftools.sage.context;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.primitives.Longs;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class RefContextCandidate implements RefContext {

    private final String sample;
    private final String chromosome;
    private final long position;
    private final Map<String, AltContextCandidate> alts;

    private int rawDepth;
    private int rawSupportRef;
    private int rawBaseQualityRef;

    public RefContextCandidate(final String sample, final String chromosome, final long position) {
        this.sample = sample;
        this.chromosome = chromosome;
        this.position = position;
        this.alts = new HashMap<>();
    }

    @NotNull
    public Collection<AltContext> alts() {
        return new ArrayList<>(alts.values());
    }

    public void refRead(int baseQuality) {
        this.rawSupportRef++;
        this.rawDepth++;
        this.rawBaseQualityRef += baseQuality;
    }

    @NotNull
    public AltContextCandidate altContext(@NotNull final String ref, @NotNull final String alt) {
        final String refAltKey = ref + "|" + alt;
        return alts.computeIfAbsent(refAltKey, key -> new AltContextCandidate(this, ref, alt));
    }


    public void altReadFixed(@NotNull final String ref, @NotNull final String alt, int baseQuality) {
        throw new UnsupportedOperationException();
    }


    public void altReadCandidate(@NotNull final String ref, @NotNull final String alt, int baseQuality, @NotNull final ReadContext interimReadContext) {
        this.rawDepth++;
        final AltContextCandidate altContext = altContext(ref, alt);
        altContext.incrementAltRead(baseQuality);
        altContext.addReadContext(interimReadContext);
    }

    @NotNull
    @Override
    public String chromosome() {
        return chromosome;
    }

    @Override
    public long position() {
        return position;
    }

    public int rawDepth() {
        return rawDepth;
    }

    @NotNull
    public String sample() {
        return sample;
    }

    public int rawSupportRef() {
        return rawSupportRef;
    }

    public int rawBaseQualityRef() {
        return rawBaseQualityRef;
    }

    @Override
    public boolean equals(@Nullable Object another) {
        if (this == another) {
            return true;
        }
        return another instanceof RefContext && equalTo((RefContext) another);
    }

    private boolean equalTo(RefContext another) {
        return chromosome().equals(another.chromosome()) && position() == another.position();
    }

    @Override
    public int hashCode() {
        int h = 5381;
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Longs.hashCode(position());
        return h;
    }
}
