package com.hartwig.hmftools.common.perf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.atomic.AtomicBoolean;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

public class StackSampler extends Thread implements AutoCloseable
{
    private final int mSleepMillis;
    private final File mOutFile;
    private final boolean mIncludeTopLineNumber;
    private final boolean mCollapseThreads;
    private final AtomicBoolean mStopSignal = new AtomicBoolean();
    private final HashMultiset<String> mStackCounts = HashMultiset.create();

    public StackSampler(final int sampleFreq, final File outFile, final boolean includeTopLineNumber, final boolean collapseThreads)
    {
        mSleepMillis = 1_000 / sampleFreq;
        mOutFile = outFile;
        mIncludeTopLineNumber = includeTopLineNumber;
        mCollapseThreads = collapseThreads;

        start();
    }

    public StackSampler(final int sampleFreq, final File outFile)
    {
        this(sampleFreq, outFile, true, true);
    }

    private void processStack(final String threadName, final StackTraceElement[] stack)
    {
        if(!threadName.toLowerCase().matches("^(main|thread-.*|gc_thread.*|g1_conc.*)$"))
            return;

        String collapsedThreadName = threadName;
        if(mCollapseThreads)
        {
            if(threadName.toLowerCase().startsWith("thread-"))
            {
                collapsedThreadName = "thread";
            }
            else if(threadName.toLowerCase().startsWith("gc_thread"))
            {
                collapsedThreadName = "gc_thread";
            }
            else if(threadName.toLowerCase().startsWith("g1_conc"))
            {
                collapsedThreadName = "g1_conc";
            }
        }

        StringJoiner foldedStack = new StringJoiner(";");
        foldedStack.add(collapsedThreadName);
        for(int i = stack.length - 1; i >= 0; i--)
        {
            StackTraceElement stackEl = stack[i];
            String stackElStr = stackEl.toString();
            String[] components;
            StringBuilder builder;
            if(i > 0 || !mIncludeTopLineNumber)
            {
                components = stackElStr.split("\\(");
                builder = new StringBuilder();
                for(int j = 0; j < components.length - 1; j++)
                    builder.append(components[j]);

                stackElStr = builder.toString();
            }

            components = stackElStr.split("/");
            builder = new StringBuilder();
            if(components.length >= 2 && !components[1].isEmpty())
                builder.append(components[1]);

            for(int j = 2; j < components.length; j++)
            {
                if(!builder.isEmpty())
                    builder.append("/");

                builder.append(components[j]);
            }

            stackElStr = builder.toString();
            if(i == 0 && stackElStr.startsWith("java.lang.Object.wait0"))
                return;

            foldedStack.add(stackElStr);
        }

        mStackCounts.add(foldedStack.toString());
    }

    @Override
    public void run()
    {
        try
        {
            final long thisThreadId = Thread.currentThread().getId();
            while(!mStopSignal.get())
            {
                Thread.sleep(mSleepMillis);
                Map<Thread, StackTraceElement[]> stackTraces = Thread.getAllStackTraces();
                for(Map.Entry<Thread, StackTraceElement[]> entry : stackTraces.entrySet())
                {
                    Thread thread = entry.getKey();
                    StackTraceElement[] stack = entry.getValue();
                    if(thisThreadId == thread.getId())
                        continue;

                    processStack(thread.getName(), stack);
                }
            }

            try(BufferedWriter writer = new BufferedWriter(new FileWriter(mOutFile)))
            {
                for(Multiset.Entry<String> entry : mStackCounts.entrySet())
                {
                    writer.write(entry.getElement() + " " + entry.getCount());
                    writer.newLine();
                }
            }
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    public void finish()
    {
        mStopSignal.set(true);
        try
        {
            join();
        }
        catch(InterruptedException e)
        {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void close() throws Exception
    {
        finish();
    }
}
