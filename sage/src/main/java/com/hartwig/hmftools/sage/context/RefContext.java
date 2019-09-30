package com.hartwig.hmftools.sage.context;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.count.ReadContext;

import org.jetbrains.annotations.NotNull;

public class RefContext implements GenomePosition {

    private final String sample;
    private final String chromosome;
    private final long position;
    private final String ref;
    private final Map<String, AltContext> alts;


    private int readDepth;
    private int mapQuality;
    private int baseQuality;
    private int quality;
    private int refReads;
    private int subprimeReadDepth;

    public RefContext(final String sample, final String chromosome, final long position, final String ref) {
        this.sample = sample;
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alts = Maps.newHashMap();
    }

    public boolean isAltsEmpty() {
        return alts.isEmpty();
    }

    public void subprimeRead(int mapQuality) {
        this.subprimeReadDepth++;
    }

    @NotNull
    public Collection<AltContext> alts() {
        return alts.values();
    }

    public void refRead(int mapQuality, int baseQuality) {
        this.readDepth++;
        this.refReads++;
        this.mapQuality += mapQuality;
        this.baseQuality += baseQuality;
        this.quality += Math.min(mapQuality, baseQuality);
    }

    @NotNull
    public AltContext altContext(@NotNull final String alt) {
        return  alts.computeIfAbsent(alt, key -> new AltContext(RefContext.this, key));
    }

    public void altRead(@NotNull final String alt, int mapQuality, int baseQuality, int recordDistance, int alignmentDistance, @NotNull final ReadContext readContext) {
        this.readDepth++;
        final AltContext altContext = alts.computeIfAbsent(alt, key -> new AltContext(RefContext.this, key));
        altContext.altRead(mapQuality, baseQuality, recordDistance, alignmentDistance,readContext);
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

    @NotNull
    public String ref() {
        return ref;
    }

    public int refQuality() {
        return quality;
    }

    public int refReads() {
        return refReads;
    }

    public int readDepth() {
        return readDepth;
    }

    public int subprimeReadDepth() {
        return subprimeReadDepth;
    }

    public String sample() {
        return sample;
    }
}
