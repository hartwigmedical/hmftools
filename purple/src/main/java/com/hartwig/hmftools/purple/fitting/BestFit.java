package com.hartwig.hmftools.purple.fitting;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BestFit
{
    public final FittedPurity Fit;
    public final FittedPurityScore Score;
    public final FittedPurityMethod Method;
    public final List<FittedPurity> AllFits;

    public BestFit(final FittedPurity fit, final FittedPurityScore score, final FittedPurityMethod method, final List<FittedPurity> allFits)
    {
        Fit = fit;
        Score = score;
        Method = method;
        AllFits = allFits;
    }

    public static List<FittedPurity> mostDiploidPerPurity(final List<FittedPurity> all)
    {
        final List<FittedPurity> sortableList = Lists.newArrayList(all);
        sortableList.sort((o1, o2) ->
        {
            double o1DistanceFromDiploid = Math.abs(2 - o1.ploidy());
            double o2DistanceFromDiploid = Math.abs(2 - o2.ploidy());
            return Double.compare(o1DistanceFromDiploid, o2DistanceFromDiploid);
        });

        final List<FittedPurity> result = Lists.newArrayList();
        final TreeSet<Double> purities = new TreeSet<>(Doubles.comparator());
        for(FittedPurity fittedPurity : sortableList)
        {
            if(purities.add(fittedPurity.purity()))
            {
                result.add(fittedPurity);
            }
        }

        Collections.sort(result);
        return result;
    }
}
