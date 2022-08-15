package com.hartwig.hmftools.common.purple;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface BestFit {

    @NotNull
    FittedPurity fit();

    @NotNull
    FittedPurityScore score();

    @NotNull
    FittedPurityMethod method();

    @NotNull
    default List<FittedPurity> bestFitPerPurity() {
        return bestFitPerPurity(allFits());
    }

    @NotNull
    List<FittedPurity> allFits();

    static List<FittedPurity> bestFitPerPurity(@NotNull final List<FittedPurity> all) {
        final List<FittedPurity> sortableList = Lists.newArrayList(all);
        sortableList.sort(Comparator.comparingDouble(FittedPurity::score));

        final List<FittedPurity> result = Lists.newArrayList();
        final TreeSet<Double> purities = new TreeSet<>(Doubles.comparator());
        for (FittedPurity fittedPurity : sortableList) {
            if (purities.add(fittedPurity.purity())) {
                result.add(fittedPurity);
            }
        }

        return result;
    }

    static List<FittedPurity> mostDiploidPerPurity(@NotNull final List<FittedPurity> all) {
        final List<FittedPurity> sortableList = Lists.newArrayList(all);
        sortableList.sort((o1, o2) -> {
            double o1DistanceFromDiploid = Math.abs(2 - o1.ploidy());
            double o2DistanceFromDiploid = Math.abs(2 - o2.ploidy());
            return Double.compare(o1DistanceFromDiploid, o2DistanceFromDiploid);
        });

        final List<FittedPurity> result = Lists.newArrayList();
        final TreeSet<Double> purities = new TreeSet<>(Doubles.comparator());
        for (FittedPurity fittedPurity : sortableList) {
            if (purities.add(fittedPurity.purity())) {
                result.add(fittedPurity);
            }
        }

        Collections.sort(result);
        return result;
    }

}
