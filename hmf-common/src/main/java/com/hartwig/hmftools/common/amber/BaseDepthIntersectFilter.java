package com.hartwig.hmftools.common.amber;

import java.util.Collection;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public class BaseDepthIntersectFilter implements Predicate<BaseDepth> {

    private boolean additional = false;
    private final Set<AmberSite> intersection = Sets.newHashSet();

    public void additional(@NotNull final Collection<BaseDepth> depth) {
        final Set<AmberSite> set = asSet(depth);
        if (additional) {
            intersection.retainAll(set);
        } else {
            additional = true;
            intersection.addAll(set);
        }
    }

    public int size() {
        return intersection.size();
    }

    @Override
    public boolean test(final BaseDepth baseDepth) {
        return !additional || intersection.contains(AmberSiteFactory.asSite(baseDepth));
    }

    @NotNull
    private static Set<AmberSite> asSet(@NotNull final Collection<BaseDepth> depth) {
        return depth.stream().map(AmberSiteFactory::asSite).collect(Collectors.toSet());
    }

}
