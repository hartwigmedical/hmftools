package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

public class ComplexCoverageCalculator
{
    private final LilacConfig mConfig;

    public ComplexCoverageCalculator(final LilacConfig config)
    {
        mConfig = config;
    }

    public List<ComplexCoverage> calculateComplexCoverages(final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes)
    {
        List<HlaAllele> alleles = Lists.newArrayList();
        complexes.stream().forEach(x -> x.Alleles.stream().filter(y -> !alleles.contains(y)).forEach(y -> alleles.add(y)));
        FragmentAlleleMatrix fragAlleleMatrix = new FragmentAlleleMatrix(fragmentAlleles, alleles);

        if(mConfig.Threads == 1 || complexes.size() < 10000) // no point in allocating to threads if complex count is small
        {
            CoverageCalcTask calcTask  = new CoverageCalcTask(0, complexes, fragAlleleMatrix, mConfig.TopScoreThreshold);
            calcTask.call();
            return calcTask.getCoverageResults();
        }

        LL_LOGGER.debug("built fragment allele matrix: fragAlleles({}) complexes({}) alleles({})",
                fragmentAlleles.size(), complexes.size(), alleles.size());

        return calcMultiThreadResults(complexes, fragAlleleMatrix);
    }

    private List<ComplexCoverage> calcMultiThreadResults(final List<HlaComplex> complexes, final FragmentAlleleMatrix fragAlleleMatrix)
    {
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("Lilac-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        List<CoverageCalcTask> coverageCalcTasks = Lists.newArrayList();
        List<FutureTask<Long>> taskList = new ArrayList<>();

        List<HlaComplex>[] complexLists = new List[mConfig.Threads];

        for(int i = 0; i < mConfig.Threads; ++i)
        {
            complexLists[i] = Lists.newArrayList();
        }

        int listIndex = 0;
        for(HlaComplex complex : complexes)
        {
            complexLists[listIndex].add(complex);

            ++listIndex;
            if(listIndex >= complexLists.length)
                listIndex = 0;
        }

        int threadIndex = 0;
        for(List<HlaComplex> complexList : complexLists)
        {
            CoverageCalcTask coverageTask = new CoverageCalcTask(threadIndex++, complexList, fragAlleleMatrix, mConfig.TopScoreThreshold);
            coverageCalcTasks.add(coverageTask);

            FutureTask<Long> futureTask = new FutureTask<>(coverageTask);

            taskList.add(futureTask);
            executorService.execute(futureTask);
        }

        try
        {
            for(FutureTask<Long> futureTask : taskList)
            {
                futureTask.get();
            }

            executorService.shutdown();

            if(mConfig.LogPerfCalcs)
            {
                PerformanceCounter combinedPerfCounter = coverageCalcTasks.get(0).getPerfCounter();

                for(int i = 1; i < coverageCalcTasks.size(); ++i)
                {
                    combinedPerfCounter.merge(coverageCalcTasks.get(i).getPerfCounter());
                }

                combinedPerfCounter.logStats();
            }

            List<ComplexCoverage> results = Lists.newArrayList();
            coverageCalcTasks.forEach(x -> results.addAll(x.getCoverageResults()));
            return results;
        }
        catch (Exception e)
        {
            LL_LOGGER.error("task execution error: {}", e.toString());
            e.printStackTrace();
            executorService.shutdown();

            return null;
        }
    }

}
