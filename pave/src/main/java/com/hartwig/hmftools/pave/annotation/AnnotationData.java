package com.hartwig.hmftools.pave.annotation;

import java.util.List;

import com.google.common.collect.Lists;

public abstract class AnnotationData
{
    public abstract String type();

    public abstract void onChromosomeComplete(final String chromosome);

    public abstract boolean enabled();

    public abstract boolean hasValidData();

    protected final List<String> mInitialChromosomes = Lists.newArrayList();

    public void registerInitialChromosomes(final List<String> chromosomes) { mInitialChromosomes.addAll(chromosomes); }
}
