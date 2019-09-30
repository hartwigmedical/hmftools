package com.hartwig.hmftools.sage.context;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.count.ReadContext;

import org.jetbrains.annotations.NotNull;

public class RefContext implements GenomePosition {

    private final String chromosome;
    private final long position;
    private final String ref;
    private final List<AltContext> alts;

    private int readDepth;
    private int refQuality;
    private int refReads;
    private int subprimeReadDepth;

    public RefContext(final String chromosome, final long position, final String ref) {
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alts = Lists.newArrayList();
    }

    public void refRead(int mapQuality, int baseQuality) {

    }

    public void altRead(@NotNull final String alt, int mapQuality, int baseQuality, @NotNull final ReadContext readContext) {

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
        return refQuality;
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
}
