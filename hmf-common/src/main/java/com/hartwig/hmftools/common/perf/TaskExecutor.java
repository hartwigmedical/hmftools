package com.hartwig.hmftools.common.perf;

import static java.lang.Thread.State.NEW;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigItemType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class TaskExecutor
{
    private TaskExecutor() {}

    static final Logger LOGGER = LogManager.getLogger(TaskExecutor.class);

    public static final String THREADS = "threads";
    public static final String THREADS_DESC = "Number of threads, 0 or 1 not multi-threaded";
    private static final int DEFAULT_THREAD_COUNT = 1;

    public static void addThreadOptions(final ConfigBuilder configBuilder)
    {
        addThreadOptions(configBuilder, DEFAULT_THREAD_COUNT);
    }

    public static void addThreadOptions(final ConfigBuilder configBuilder, final int defaultCount)
    {
        configBuilder.addConfigItem(ConfigItemType.INTEGER, THREADS, false, THREADS_DESC, String.valueOf(defaultCount));

        setDefaultThreadExceptionHandler();
    }

    public static int parseThreads(final ConfigBuilder configBuilder) { return configBuilder.getInteger(THREADS); }

    public static <T> boolean executeTasks(final List<Callable<T>> tasks, int threadCount)
    {
        if(threadCount <= 1)
        {
            for(Callable<T> task : tasks)
            {
                try
                {
                    task.call();
                }
                catch(Exception e)
                {
                    LOGGER.error("task execution error: {}", e.toString());
                    e.printStackTrace();
                    return false;
                }
            }

            return true;
        }

        int digits = Integer.toString(threadCount - 1).length();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + digits + "d").build();

        ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);
        List<Future<T>> threadTaskList = new ArrayList<>();

        for(Callable<T> task : tasks)
        {
            threadTaskList.add(executorService.submit(task));
        }

        if(!checkTaskCompletion(threadTaskList))
        {
            LOGGER.info("shutting down remaining tasks");
            threadTaskList.forEach(x -> x.cancel(true));
            executorService.shutdown();
            return false;
        }

        executorService.shutdown();
        return true;
    }

    public static boolean executeRunnables(final List<? extends Runnable> tasks, int threadCount)
    {
        return executeTasks(tasks.stream().map(Executors::callable).collect(Collectors.toList()), threadCount);
    }

    private static <T> boolean checkTaskCompletion(final List<Future<T>> taskList)
    {
        try
        {
            for(Future<T> future : taskList)
            {
                future.get();
            }
        }
        catch(Exception e)
        {
            LOGGER.error("task execution error: {}", e.toString());
            e.printStackTrace();
            return false;
        }

        return true;
    }

    public static boolean runThreadTasks(final List<? extends Thread> workers)
    {
        workers.stream().filter(x -> x.getState() == NEW).forEach(Thread::start);

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                LOGGER.error("thread execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        return true;
    }

    public static void setDefaultThreadExceptionHandler()
    {
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            LOGGER.error("thread exception: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        });
    }
}
