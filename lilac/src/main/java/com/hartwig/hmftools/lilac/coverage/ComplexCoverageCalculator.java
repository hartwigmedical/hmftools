package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

public class ComplexCoverageCalculator
{
    private final int mThreadCount;
    private final double mTopScoreThreshold;

    public ComplexCoverageCalculator(int threadCount, double topScoreThreshold)
    {
        mThreadCount = threadCount;
        mTopScoreThreshold = topScoreThreshold;
    }

    public List<ComplexCoverage> calculateComplexCoverages(final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes)
    {
        List<HlaAllele> alleles = Lists.newArrayList();
        complexes.stream().forEach(x -> x.Alleles.stream().filter(y -> !alleles.contains(y)).forEach(y -> alleles.add(y)));
        FragmentAlleleMatrix fragAlleleMatrix = new FragmentAlleleMatrix(fragmentAlleles, alleles);

        if(mThreadCount > 1 && complexes.size() >= 10000) // no point in allocating to threads if complex count is small
        {
            return calcMultiThreadResults(complexes, fragAlleleMatrix);
        }

        CoverageCalcTask calcTask  = new CoverageCalcTask(0, complexes, fragAlleleMatrix, mTopScoreThreshold);
        calcTask.call();
        return calcTask.getCoverageResults();
    }

    private List<ComplexCoverage> calcMultiThreadResults(final List<HlaComplex> complexes, final FragmentAlleleMatrix fragAlleleMatrix)
    {
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("Lilac-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mThreadCount, namedThreadFactory);

        List<CoverageCalcTask> coverageCalcTasks = Lists.newArrayList();
        List<FutureTask<Long>> taskList = new ArrayList<>();

        List<HlaComplex>[] complexLists = new List[mThreadCount];

        for(int i = 0; i < mThreadCount; ++i)
        {
            complexLists[i] = Lists.newArrayList();
        }

        int listIndex = 0;
        for (HlaComplex complex : complexes)
        {
            complexLists[listIndex].add(complex);

            ++listIndex;
            if(listIndex >= complexLists.length)
                listIndex = 0;
        }

        int threadIndex = 0;
        for(List<HlaComplex> complexList : complexLists)
        {
            CoverageCalcTask coverageTask = new CoverageCalcTask(threadIndex++, complexList, fragAlleleMatrix, mTopScoreThreshold);
            coverageCalcTasks.add(coverageTask);

            FutureTask<Long> futureTask = new FutureTask<>(coverageTask);

            taskList.add(futureTask);
            executorService.execute(futureTask);
        }

        try
        {
            for (FutureTask<Long> futureTask : taskList)
            {
                futureTask.get();
            }

            executorService.shutdown();

            coverageCalcTasks.forEach(x -> x.logPerfResults());

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
