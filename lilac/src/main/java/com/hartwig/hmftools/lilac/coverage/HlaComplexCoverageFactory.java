package com.hartwig.hmftools.lilac.coverage;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.progress.FutureProgressTracker;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import org.jetbrains.annotations.NotNull;

public class HlaComplexCoverageFactory
{
    private final FutureProgressTracker mProgressTracker;
    private final LilacConfig mConfig;

    public HlaComplexCoverageFactory(@NotNull LilacConfig config)
    {
        mConfig = config;
        mProgressTracker = new FutureProgressTracker(0.1, 10000);
    }

    public final List<HlaAllele> rankedGroupCoverage(
            int take, @NotNull List<FragmentAlleles> fragmentAlleles,
            @NotNull List<HlaComplex> complexes)
    {
        return Lists.newArrayList();
        
        /*
        void $receiver$iv$iv;
        Comparable<HlaComplexCoverage> comparable;
        Comparable<HlaComplexCoverage> it;
        Object object;
        Iterable $receiver$iv$iv2;
        Iterable $receiver$iv;
        Intrinsics.checkParameterIsNotNull(fragmentAlleles, (String) "fragmentAlleles");
        Intrinsics.checkParameterIsNotNull(complexes, (String) "complexes");
        Iterable iterable = complexes;
        void var6_5 = $receiver$iv;
        Object destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv2)
        {
            HlaComplex hlaComplex = (HlaComplex) item$iv$iv;
            object = destination$iv$iv;
            boolean bl = false;
            comparable = Companion.proteinCoverage(fragmentAlleles, (Collection<HlaAllele>) ((HlaComplex) ((Object) it)).getAlleles());
            object.add(comparable);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv2 = $receiver$iv;
        destination$iv$iv = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                HlaComplexCoverage it = (HlaComplexCoverage) a;
                boolean bl = false;
                Comparable comparable = Integer.valueOf(-it.getTotalCoverage());
                it = (HlaComplexCoverage) b;
                Comparable comparable2 = comparable;
                bl = false;
                Integer n = -it.getTotalCoverage();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) n);
            }
        };
        $receiver$iv = CollectionsKt.sortedWith((Iterable) $receiver$iv$iv2, (Comparator) destination$iv$iv);
        $receiver$iv$iv2 = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            it = (HlaComplexCoverage) element$iv$iv;
            boolean bl = false;
            Iterable list$iv$iv = ((HlaComplexCoverage) it).getAlleleCoverage();
            CollectionsKt.addAll((Collection) destination$iv$iv, (Iterable) list$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv2 = $receiver$iv;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv2)
        {
            it = (HlaAlleleCoverage) item$iv$iv;
            object = destination$iv$iv;
            boolean bl = false;
            comparable = ((HlaAlleleCoverage) it).getAllele();
            object.add(comparable);
        }
        List topRanked = CollectionsKt.distinct((Iterable) ((List) destination$iv$iv));
        List topTakers = CollectionsKt.take((Iterable) topRanked, (int) take);
        Iterable $receiver$iv2 = topRanked;
        Iterable $i$f$sortedBy = $receiver$iv2;
        Collection destination$iv$iv2 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            HlaAllele it2 = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!this.mConfig.getCommonAlleles().contains(it2))
            {
                continue;
            }
            destination$iv$iv2.add(element$iv$iv);
        }
        List topRankedKeepers = (List) destination$iv$iv2;
        return CollectionsKt.distinct((Iterable) CollectionsKt.plus((Collection) topTakers, (Iterable) topRankedKeepers));
        
         */
    }

    public final List<HlaComplexCoverage> rankedComplexCoverage(@NotNull ExecutorService executorService,
            @NotNull List<FragmentAlleles> fragmentAlleles, @NotNull List<HlaComplex> complexes, @NotNull List<HlaAllele> recovered)
    {
        return Lists.newArrayList();

        /*
        void $receiver$iv$iv;
        Object untrackedCallable2;
        Intrinsics.checkParameterIsNotNull((Object) executorService, (String) "executorService");
        Intrinsics.checkParameterIsNotNull(fragmentAlleles, (String) "fragmentAlleles");
        Intrinsics.checkParameterIsNotNull(complexes, (String) "complexes");
        Intrinsics.checkParameterIsNotNull(recovered, (String) "recovered");
        HlaComplexCoverageRanking ranking =
                new HlaComplexCoverageRanking(this.mConfig.getMaxDistanceFromTopScore(), this.mConfig.getCommonAlleles(), recovered, this.mConfig
                        .getStopLossRecoveryAlleles());
        List list = new ArrayList();
        for(HlaComplex complex : complexes)
        {
            Callable trackedCallable;
            untrackedCallable2 = new Callable<HlaComplexCoverage>(fragmentAlleles, complex)
            {
                public final HlaComplexCoverage call()
                {
                    return HlaComplexCoverageFactory.Companion.proteinCoverage(this.$fragmentAlleles, (Collection<HlaAllele>) this.$complex.getAlleles());
                }

                {
                    this.$fragmentAlleles = list;
                    this.$complex = hlaComplex;
                }
            };
            Intrinsics.checkExpressionValueIsNotNull((Object) this.mProgressTracker.add((Callable) untrackedCallable2), (String) "progressTracker.add(untrackedCallable)");
            Future future = executorService.submit(trackedCallable);
            Intrinsics.checkExpressionValueIsNotNull(future, (String) "executorService.submit(trackedCallable)");
            list.add(future);
        }
        Iterable $receiver$iv = list;
        untrackedCallable2 = $receiver$iv;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            Future future = (Future) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            HlaComplexCoverage
                    hlaComplexCoverage = (HlaComplexCoverage) it.get();
            collection.add(hlaComplexCoverage);
        }
        List result = (List) destination$iv$iv;
        return ranking.candidateRanking(result);

         */
    }

    public final HlaComplexCoverage groupCoverage(@NotNull List<FragmentAlleles> fragmentAlleles,
            @NotNull Collection<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = this.fragmentAlleles(fragmentAlleles, alleles);
        return null;
        //return HlaComplexCoverage.create(HlaAlleleCoverage.groupCoverage(filteredFragments));
    }

    @NotNull
    public final HlaComplexCoverage proteinCoverage(
            @NotNull List<FragmentAlleles> fragmentAlleles,
            @NotNull Collection<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = this.fragmentAlleles(fragmentAlleles, alleles);
        return null;
        //return HlaComplexCoverage.create(HlaAlleleCoverage.proteinCoverage(filteredFragments));
    }

    private final List<FragmentAlleles> fragmentAlleles(List<FragmentAlleles> fragmentAlleles, Collection<HlaAllele> alleles)
    {
        return Lists.newArrayList();
        //return FragmentAlleles.Companion.filter(fragmentAlleles, alleles);
    }
}
