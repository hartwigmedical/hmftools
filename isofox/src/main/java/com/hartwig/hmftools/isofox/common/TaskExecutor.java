package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

public class TaskExecutor
{
    public static boolean executeChromosomeTask(final List<Callable> tasks, int threadCount)
    {
        // chrTasks.forEach(x -> x.setTaskType(taskType));

        if(threadCount <= 1)
        {
            for(Callable task : tasks)
            {
                try
                {
                    task.call();
                }
                catch(Exception e)
                {
                    ISF_LOGGER.error("task execution error: {}", e.toString());
                }
            }

            return true;
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("Isofox-%d").build();

        ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);
        List<FutureTask> threadTaskList = new ArrayList<FutureTask>();

        for(Callable task : tasks)
        {
            FutureTask futureTask = new FutureTask(task);

            threadTaskList.add(futureTask);
            executorService.execute(futureTask);
        }

        if(!checkThreadCompletion(threadTaskList))
            return false;

        executorService.shutdown();
        return true;
    }

    private static boolean checkThreadCompletion(final List<FutureTask> taskList)
    {
        try
        {
            for (FutureTask futureTask : taskList)
            {
                futureTask.get();
            }
        }
        catch (Exception e)
        {
            ISF_LOGGER.error("task execution error: {}", e.toString());
            e.printStackTrace();
            return false;
        }

        return true;
    }

}
