package com.hartwig.hmftools.orange.algo.wildtype;

import com.hartwig.hmftools.datamodel.wildtype.ImmutableWildTypeGene;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;

import org.jetbrains.annotations.NotNull;

public final class TestWildTypeFactory
{

    private TestWildTypeFactory()
    {
    }

    @NotNull
    public static WildTypeGene create(@NotNull String gene)
    {
        return ImmutableWildTypeGene.builder().gene(gene).build();
    }
}
