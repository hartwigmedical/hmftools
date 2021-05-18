package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.progress.FutureProgressTracker;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

public class ComplexCoverageCalculator
{
    private final FutureProgressTracker mProgressTracker;
    private final LilacConfig mConfig;

    public ComplexCoverageCalculator(final LilacConfig config)
    {
        mConfig = config;
        mProgressTracker = new FutureProgressTracker(0.1, 10000);
    }

    public List<HlaComplexCoverage> calculateComplexCoverages(final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes)
    {
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("Lilac-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

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

            executorService.shutdown();

            List<HlaComplexCoverage> results = coverageCalcTasks.stream().map(x -> x.getCoverage()).collect(Collectors.toList());
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
