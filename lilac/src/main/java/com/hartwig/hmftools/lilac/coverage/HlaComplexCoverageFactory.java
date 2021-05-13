package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.coverage.CoverageCalcTask.proteinCoverage;
import static com.hartwig.hmftools.lilac.hla.HlaAllele.contains;
import static com.hartwig.hmftools.lilac.hla.HlaAllele.dedup;
import static com.hartwig.hmftools.lilac.hla.HlaAllele.takeN;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.progress.FutureProgressTracker;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.FutureTask;
import java.util.stream.Collectors;

public class HlaComplexCoverageFactory
{
    private final FutureProgressTracker mProgressTracker;
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    public HlaComplexCoverageFactory(final LilacConfig config, final ReferenceData refData)
    {
        mConfig = config;
        mRefData = refData;
        mProgressTracker = new FutureProgressTracker(0.1, 10000);
    }

    public List<HlaAllele> rankedGroupCoverage(
            int take, final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes)
    {
        List<HlaComplexCoverage> coverages = complexes.stream()
                .map(x -> proteinCoverage(fragmentAlleles, x.getAlleles())).collect(Collectors.toList());

        Collections.sort(coverages, new HlaComplexCoverage.TotalCoverageSorter());
        List<HlaAllele> topRanked = Lists.newArrayList();

        for(HlaComplexCoverage coverage : coverages)
        {
            coverage.getAlleleCoverage().stream().filter(x -> !topRanked.contains(x.Allele)).forEach(x -> topRanked.add(x.Allele));
        }

        List<HlaAllele> topTakers = takeN(topRanked, take);

        List<HlaAllele> topRankedKeepers = topRanked.stream().filter(x -> contains(mRefData.CommonAlleles, x)).collect(Collectors.toList());

        List<HlaAllele> uniqueRanked = topTakers.stream().collect(Collectors.toList());
        topRankedKeepers.stream().filter(x -> !contains(topTakers, x)).forEach(x -> uniqueRanked.add(x));

        return uniqueRanked;
    }

    public List<HlaComplexCoverage> rankedComplexCoverage(
            final ExecutorService executorService,
            final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes, final List<HlaAllele> recovered)
    {
        HlaComplexCoverageRanking ranking = new HlaComplexCoverageRanking(
                mConfig.MaxDistanceFromTopScore, mRefData.CommonAlleles, recovered, mRefData.StopLossRecoveryAlleles);

        List<CoverageCalcTask> coverageCalcTasks = Lists.newArrayList();
        List<FutureTask> taskList = new ArrayList<FutureTask>();

        for (HlaComplex complex : complexes)
        {
            CoverageCalcTask coverageTask = new CoverageCalcTask(fragmentAlleles, complex);
            coverageCalcTasks.add(coverageTask);

            FutureTask futureTask = new FutureTask(coverageTask);

            taskList.add(futureTask);
            executorService.execute(futureTask);

            // mProgressTracker.add(futureTask);
            // taskList.add(executorService.submit(coverageTask));

            /*
            val untrackedCallable: Callable<HlaComplexCoverage> = Callable { proteinCoverage(fragmentAlleles, complex.alleles) }
            val trackedCallable: Callable<HlaComplexCoverage> = progressTracker.add(untrackedCallable)
            list.add(executorService.submit(trackedCallable))
             */
        }

        try
        {
            for (FutureTask futureTask : taskList)
            {
                futureTask.get();
            }

            List<HlaComplexCoverage> results = coverageCalcTasks.stream().map(x -> x.getCoverage()).collect(Collectors.toList());
            return ranking.candidateRanking(results);

        }
        catch (Exception e)
        {
            LL_LOGGER.error("task execution error: {}", e.toString());
            e.printStackTrace();
            return null;
        }
    }
}
