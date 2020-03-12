package com.hartwig.hmftools.sage.context;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.primitives.Longs;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class RefContextFixed implements RefContext {

    private final String sample;
    private final String chromosome;
    private final long position;
    private final Map<String, AltContextFixed> alts;

    private int rawDepth;
    private int rawSupportRef;
    private int rawBaseQualityRef;

    RefContextFixed(final String sample, final String chromosome, final long position) {
        this.sample = sample;
        this.chromosome = chromosome;
        this.position = position;
        this.alts = new HashMap<>();
    }

    @NotNull
    public Collection<AltContext> alts() {
        return new ArrayList<>(alts.values());
    }

    @NotNull
    public Collection<AltContextFixed> fixedAlts() {
        return alts.values();
    }

    public void refRead(int baseQuality) {
        this.rawSupportRef++;
        this.rawDepth++;
        this.rawBaseQualityRef += baseQuality;
    }

    @Nullable
    public AltContextFixed altContext(@NotNull final String ref, @NotNull final String alt) {
        final String refAltKey = ref + "|" + alt;
        return alts.get(refAltKey);
    }

    public AltContext altContext(@NotNull final String ref, @NotNull final String alt, boolean realign, @NotNull final ReadContext readContext) {
        final String refAltKey = ref + "|" + alt;
        return alts.computeIfAbsent(refAltKey, key -> new AltContextFixed(this, ref, alt, realign, readContext));
    }

    public void altReadFixed(@NotNull final String ref, @NotNull final String alt, int baseQuality) {
        this.rawDepth++;
        final AltContextFixed altContext = altContext(ref, alt);
        if (altContext != null) {
            altContext.incrementAltRead(baseQuality);
        }
    }

    public void altReadCandidate(@NotNull final String ref, @NotNull final String alt, int baseQuality,
            @NotNull final ReadContext interimReadContext) {
        throw new UnsupportedOperationException();
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
