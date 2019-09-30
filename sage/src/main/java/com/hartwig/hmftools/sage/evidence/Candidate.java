package com.hartwig.hmftools.sage.evidence;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.count.ReadContext;
import com.hartwig.hmftools.sage.count.ReadContextConsumer;

import org.jetbrains.annotations.NotNull;

public class Candidate implements VariantHotspot {

    private final String chromosome;
    private final long position;
    private final String ref;
    private final String alt;

    private int readDepth;
    private int refSupport;
    private int refQuality;
    private int altSupport;
    private int altQuality;

    private final Map<ReadContext, ReadContextConsumer> interimCounts;

    public Candidate(final String chromosome, final long position, final String ref, final String alt) {
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
        this.interimCounts = Maps.newHashMap();
    }

    public void addSupport(int quality, @NotNull final ReadContext context) {
        refSupport += 1;
        refQuality += quality;

        if (!interimCounts.keySet().contains(context)) {
            interimCounts.put(context, new ReadContextConsumer(this, context));
        }

        interimCounts.values().forEach(x -> x.accept(position, context));
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
    @Override
    public String ref() {
        return ref;
    }

    @NotNull
    @Override
    public String alt() {
        return alt;
    }
}
