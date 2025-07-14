package com.hartwig.hmftools.common.perf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.time.Duration;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.pcollections.HashTreePMap;
import org.pcollections.PMap;

public class StackSampler implements AutoCloseable
{
    private static final Logger LOGGER = LogManager.getLogger(StackSampler.class);

    public static final String STACK_SAMPLER_THREAD_NAME = "STACK_SAMPLER";
    public static final String STACK_PROCESSOR_THREAD_NAME = "STACK_PROCESSOR";

    private final Duration mSamplePeriod;
    private final File mOutFile;
    private final boolean mIncludeTopLineNumber;
    private final boolean mCollapseThreads;
    private final AtomicBoolean mStopSignal = new AtomicBoolean();
    private final ConcurrentLinkedQueue<Map<Thread, StackTraceElement[]>> mSampleQueue = new ConcurrentLinkedQueue<>();
    private final AtomicReference<PMap<String, Long>> mStackCounts = new AtomicReference<>(HashTreePMap.empty());

    private Thread mStackSampler;
    private Thread mStackProcessor;

    public StackSampler(final int samplesPerSecond, final File outFile, final boolean includeTopLineNumber, final boolean collapseThreads)
    {
        mSamplePeriod = Duration.ofSeconds(1).dividedBy(samplesPerSecond);
        mOutFile = outFile;
        mIncludeTopLineNumber = includeTopLineNumber;
        mCollapseThreads = collapseThreads;

        run();
    }

    public StackSampler(final int sampleFreq, final File outFile)
    {
        this(sampleFreq, outFile, true, true);
    }

    private boolean threadNameFilter(final String threadName)
    {
        if(threadName.equals(STACK_SAMPLER_THREAD_NAME))
            return false;

        if(threadName.equals(STACK_PROCESSOR_THREAD_NAME))
            return true;

        return threadName.toLowerCase().matches("^(main|thread-.*|gc_thread.*|g1_conc.*)$");
    }

    private static String collapseThreadName(final String threadName)
    {
        if(threadName.toLowerCase().startsWith("thread-"))
        {
            return "thread";
        }

        if(threadName.toLowerCase().startsWith("gc_thread"))
        {
            return "gc_thread";
        }

        if(threadName.toLowerCase().startsWith("g1_conc"))
        {
            return "g1_conc";
        }

        return threadName;
    }

    private static String stripLineNumberFromStackEl(final String stackEl)
    {
        String[] components = stackEl.split("\\(");
        StringBuilder builder = new StringBuilder();
        for(int j = 0; j < components.length - 1; j++)
            builder.append(components[j]);

        return builder.toString();
    }

    private static String stripPrefixFromStackEl(final String stackEl)
    {
        String[] components = stackEl.split("/");
        StringBuilder builder = new StringBuilder();
        if(components.length >= 2 && !components[1].isEmpty())
            builder.append(components[1]);

        for(int j = 2; j < components.length; j++)
        {
            if(!builder.isEmpty())
                builder.append("/");

            builder.append(components[j]);
        }

        return builder.toString();
    }

    private void processStack(final String threadName_, final StackTraceElement[] stack)
    {
        if(stack.length == 0)
            return;

        if(!threadNameFilter(threadName_))
            return;

        String threadName = mCollapseThreads ? collapseThreadName(threadName_) : threadName_;

        StringJoiner foldedStack = new StringJoiner(";");
        foldedStack.add(threadName);

        List<String> stackEls = Lists.newArrayList(Arrays.stream(stack).map(StackTraceElement::toString).toList());
        Collections.reverse(stackEls);
        for(String stackEl : stackEls)
        {
            stackEl = stripPrefixFromStackEl(stackEl);
            stackEl = stripLineNumberFromStackEl(stackEl);
            foldedStack.add(stackEl);
        }

        if(mIncludeTopLineNumber)
            foldedStack.add(stack[0].toString());

        final String foldedStackStr = foldedStack.toString().replace(' ', '_');

        mStackCounts.updateAndGet(stackCounts ->
        {
            if(!stackCounts.containsKey(foldedStackStr))
                return stackCounts.plus(foldedStackStr, 1L);

            long currentCount = stackCounts.get(foldedStackStr);
            return stackCounts.plus(foldedStackStr, currentCount + 1L);
        });
    }

    private void run()
    {
        mStackSampler = new Thread(() ->
        {
            long nextSampleTime = System.nanoTime() + mSamplePeriod.toNanos();
            while(!mStopSignal.get())
            {
                long currentTime = System.nanoTime();
                if(currentTime < nextSampleTime)
                {
                    Duration sleepDuration = Duration.ofNanos(nextSampleTime - currentTime);
                    try
                    {
                        Thread.sleep(sleepDuration.toMillis());
                    }
                    catch(InterruptedException e)
                    {
                    }
                }

                nextSampleTime += mSamplePeriod.toNanos();
                mSampleQueue.add(Thread.getAllStackTraces());
            }
        }, STACK_SAMPLER_THREAD_NAME);

        mStackSampler.start();

        mStackProcessor = new Thread(() ->
        {
            while(!(mStopSignal.get() && mSampleQueue.isEmpty()))
            {
                Map<Thread, StackTraceElement[]> sample = mSampleQueue.poll();
                if(sample == null)
                {
                    try
                    {
                        Thread.sleep(1_000);
                    }
                    catch(InterruptedException e)
                    {
                    }

                    continue;
                }

                for(Map.Entry<Thread, StackTraceElement[]> entry : sample.entrySet())
                {
                    String threadName = entry.getKey().getName();
                    StackTraceElement[] stackTrace = entry.getValue();
                    processStack(threadName, stackTrace);
                }
            }
        }, STACK_PROCESSOR_THREAD_NAME);

        mStackProcessor.start();
    }

    @Override
    public void close() throws Exception
    {
        mStopSignal.set(true);
        mStackSampler.join();
        mStackProcessor.join();

        long totalStackCount = 0;
        try(BufferedWriter writer = new BufferedWriter(new FileWriter(mOutFile)))
        {
            for(Map.Entry<String, Long> entry : mStackCounts.get().entrySet())
            {
                totalStackCount += entry.getValue();
                writer.write(entry.getKey() + " " + entry.getValue());
                writer.newLine();
            }
        }

        LOGGER.info("StackSampler: {} stacks written", totalStackCount);
    }
}
