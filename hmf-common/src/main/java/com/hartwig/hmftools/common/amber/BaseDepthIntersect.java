package com.hartwig.hmftools.common.amber;

import java.util.Collection;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.jetbrains.annotations.NotNull;

public class BaseDepthIntersect {

    private final Set<AmberSite> intersection = Sets.newHashSet();

    public void add(@NotNull final Collection<BaseDepth> depth) {
        Set<AmberSite> set = asSet(depth);
        if (intersection.isEmpty()) {
            intersection.addAll(set);
        } else {
            intersection.retainAll(set);
        }
    }

    @NotNull
    public ListMultimap<Chromosome, BaseDepth> intersect(@NotNull final ListMultimap<Chromosome, BaseDepth> depth) {
        final ListMultimap<Chromosome, BaseDepth> result = ArrayListMultimap.create();
        for (Chromosome chromosome : depth.keySet()) {
            for (BaseDepth baseDepth : depth.get(chromosome)) {
                if (intersection.contains(asSite(baseDepth))) {
                    result.put(chromosome, baseDepth);
                }
            }
        }

        return result;
    }

    @NotNull
    private static Set<AmberSite> asSet(@NotNull final Collection<BaseDepth> depth) {
        return depth.stream().map(BaseDepthIntersect::asSite).collect(Collectors.toSet());
    }

    @NotNull
    private static AmberSite asSite(@NotNull final BaseDepth baseDepth) {
        return ImmutableAmberSite.builder()
                .from(baseDepth)
                .snpCheck(false)
                .ref(baseDepth.ref().toString())
                .alt(baseDepth.alt().toString())
                .build();
    }

}
