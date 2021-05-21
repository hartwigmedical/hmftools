package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.progress.FutureProgressTracker;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

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

        List<HlaComplex>[] complexLists = new List[mConfig.Threads];

        for(int i = 0; i < mConfig.Threads; ++i)
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

        Map<HlaAllele,List<FragmentAlleles>> alleleFragmentMap = null; // buildAlleleFragmentMap(fragmentAlleles, complexes);

        for(List<HlaComplex> complexList : complexLists)
        {
            CoverageCalcTask coverageTask = new CoverageCalcTask(fragmentAlleles, complexList, alleleFragmentMap);
            coverageCalcTasks.add(coverageTask);

            FutureTask futureTask = new FutureTask(coverageTask);

            taskList.add(futureTask);
            executorService.execute(futureTask);
        }

        /*
        for (HlaComplex complex : complexes)
        {
            CoverageCalcTask coverageTask = new CoverageCalcTask(fragmentAlleles, complex);
            coverageCalcTasks.add(coverageTask);

            FutureTask futureTask = new FutureTask(coverageTask);

            taskList.add(futureTask);
            executorService.execute(futureTask);

            // mProgressTracker.add(futureTask);
            // taskList.add(executorService.submit(coverageTask));


           // progress tracker usage
            // val untrackedCallable: Callable<HlaComplexCoverage> = Callable { proteinCoverage(fragmentAlleles, complex.alleles) }
            // val trackedCallable: Callable<HlaComplexCoverage> = progressTracker.add(untrackedCallable)
            // list.add(executorService.submit(trackedCallable))
        }
        */

        try
        {
            for (FutureTask futureTask : taskList)
            {
                futureTask.get();
            }

            executorService.shutdown();

            coverageCalcTasks.forEach(x -> x.logPerfResults());

            List<HlaComplexCoverage> results = Lists.newArrayList();
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

    private Map<HlaAllele,List<FragmentAlleles>> buildAlleleFragmentMap(
            final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes)
    {
        int allelesWithFrags = 0;
        Map<HlaAllele,List<FragmentAlleles>> alleleFragmentMap = Maps.newHashMap();

        for(HlaComplex complex : complexes)
        {
            for(HlaAllele allele : complex.getAlleles())
            {
                if(alleleFragmentMap.containsKey(allele))
                    continue;

                List<FragmentAlleles> fragments = Lists.newArrayList();
                for(FragmentAlleles fragment : fragmentAlleles)
                {
                    if(fragment.getFull().contains(allele) || fragment.getPartial().contains(allele))
                    {
                        fragments.add(fragment);
                    }
                }

                if(!fragments.isEmpty())
                    ++allelesWithFrags;

                alleleFragmentMap.put(allele, fragments);
            }
        }

        LL_LOGGER.info("alleles {} with fragments: {}", alleleFragmentMap.size(), allelesWithFrags);
        return alleleFragmentMap;
    }

    private static FragmentAlleles createIfContains(final FragmentAlleles fragAlleles, final HlaAllele allele)
    {
        List<HlaAllele> emptyList = Lists.newArrayList();
        List<HlaAllele> fullList = fragAlleles.getFull().contains(allele) ? Lists.newArrayList(allele) : emptyList;
        List<HlaAllele> partialList = fragAlleles.getPartial().contains(allele) ? Lists.newArrayList(allele) : emptyList;

        if(fullList.isEmpty() && partialList.isEmpty())
            return null;

        return new FragmentAlleles(fragAlleles.getFragment(), fullList, partialList, emptyList);
    }

}
