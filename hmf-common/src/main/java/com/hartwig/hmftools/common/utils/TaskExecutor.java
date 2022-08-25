package com.hartwig.hmftools.common.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TaskExecutor
{
    private static final Logger LOGGER = LogManager.getLogger(TaskExecutor.class);

    public static final String THREADS = "threads";
    public static final String THREADS_DESC = "Number of threads";

    public static void addThreadOptions(final Options options)
    {
        options.addOption(THREADS, true, THREADS_DESC);
    }

    public static int parseThreads(final CommandLine cmd) { return Integer.parseInt(cmd.getOptionValue(THREADS, "1")); }

    public static boolean executeTasks(final List<Callable> tasks, int threadCount)
    {
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
                    LOGGER.error("task execution error: {} stack: {}", e.toString());
                    e.printStackTrace();
                }
            }

            return true;
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("Thread-%d").build();

        ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);
        List<FutureTask> threadTaskList = new ArrayList<FutureTask>();

        for(Callable task : tasks)
        {
            FutureTask futureTask = new FutureTask(task);

            threadTaskList.add(futureTask);
            executorService.execute(futureTask);
        }

        if(!checkThreadCompletion(threadTaskList))
        {
            LOGGER.info("shutting down remaining tasks");
            threadTaskList.forEach(x -> x.cancel(true));
            executorService.shutdown();
            return false;
        }

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
            LOGGER.error("task execution error: {}", e.toString());
            e.printStackTrace();
            return false;
        }

        return true;
    }

}
