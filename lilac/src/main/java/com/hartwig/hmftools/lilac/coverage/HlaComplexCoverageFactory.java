package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.hla.HlaAllele.takeN;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.progress.FutureProgressTracker;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class HlaComplexCoverageFactory
{
    private final FutureProgressTracker mProgressTracker;
    private final LilacConfig mConfig;

    public HlaComplexCoverageFactory(final LilacConfig config)
    {
        mConfig = config;
        mProgressTracker = new FutureProgressTracker(0.1, 10000);
    }

    public List<HlaAllele> rankedGroupCoverage(
            int take, final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes)
    {
        List<HlaComplexCoverage> coverage = complexes.stream()
                .map(x -> proteinCoverage(fragmentAlleles, x.getAlleles())).collect(Collectors.toList());

        Collections.sort(coverage, new HlaComplexCoverage.TotalCoverageSorter());
        List<HlaAllele> coverageAlleles = Lists.newArrayList();
        coverage.forEach(x -> x.getAlleleCoverage().forEach(y -> coverageAlleles.add(y.Allele)));

        List<HlaAllele> topRanked = coverageAlleles.stream().distinct().collect(Collectors.toList());

        List<HlaAllele> topTakers = takeN(topRanked, take);

        List<HlaAllele> topRankedKeepers = topRanked.stream().filter(x -> mConfig.CommonAlleles.contains(x)).collect(Collectors.toList());

        List<HlaAllele> uniqueRanked = topTakers.stream().collect(Collectors.toList());
        topRankedKeepers.stream().filter(x -> !topTakers.contains(x)).forEach(x -> uniqueRanked.add(x));

        return uniqueRanked;
    }

    public List<HlaComplexCoverage> rankedComplexCoverage(
            final ExecutorService executorService,
            final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes, final List<HlaAllele> recovered)
    {
        HlaComplexCoverageRanking ranking = new HlaComplexCoverageRanking(
                mConfig.MaxDistanceFromTopScore, mConfig.CommonAlleles, recovered, mConfig.StopLossRecoveryAlleles);

        List<Future<HlaComplexCoverage>> futures = Lists.newArrayList();

        for (HlaComplex complex : complexes)
        {
            HlaComplexCoverage coverage = proteinCoverage(fragmentAlleles, complex.getAlleles());

            //Callable<HlaComplexCoverage> untrackedCallable = coverage;
            Callable untrackedCallable = null; //new Callable<HlaComplexCoverage>(fragmentAlleles, complex);

            Callable<HlaComplexCoverage> trackedCallable = mProgressTracker.add(untrackedCallable);
            futures.add(executorService.submit(trackedCallable));

            /*
            val untrackedCallable: Callable<HlaComplexCoverage> = Callable { proteinCoverage(fragmentAlleles, complex.alleles) }
            val trackedCallable: Callable<HlaComplexCoverage> = progressTracker.add(untrackedCallable)
            list.add(executorService.submit(trackedCallable))
             */
        }

        try
        {
            List<HlaComplexCoverage> results = Lists.newArrayList(); // )futures.stream().map(x -> x.get()).collect(Collectors.toList());
            return ranking.candidateRanking(results);
        }
        catch(Exception e)
        {
            LL_LOGGER.error("execution exception: {}", e.toString());
            return null;
        }

        /*

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

    public static HlaComplexCoverage groupCoverage(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = fragmentAlleles(fragmentAlleles, alleles);
        return HlaComplexCoverage.create(HlaAlleleCoverage.groupCoverage(filteredFragments));
    }

    public static HlaComplexCoverage proteinCoverage(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = fragmentAlleles(fragmentAlleles, alleles);
        return HlaComplexCoverage.create(HlaAlleleCoverage.proteinCoverage(filteredFragments));
    }

    private static List<FragmentAlleles> fragmentAlleles(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        return FragmentAlleles.filter(fragmentAlleles, alleles);
    }
}
