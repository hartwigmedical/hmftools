package com.hartwig.hmftools.bamtools.bamtofastq;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.bamtofastq.BamToFastqConfig.BFQ_LOGGER;
import static com.hartwig.hmftools.bamtools.bamtofastq.BamToFastqConfig.addConfig;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;

import com.google.common.io.Closer;
import com.google.inject.Guice;
import com.google.inject.Inject;
import com.google.inject.Injector;
import com.google.inject.Stage;
import com.hartwig.hmftools.bamtools.bamtofastq.readcache.ReadCacheInterface;
import com.hartwig.hmftools.bamtools.bamtofastq.writer.PairedFastqWriterInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// TODO NEXT: TEST
// TODO NEXT: By read group
// TODO: Flame graph.
// TODO: move classes to hmf common
// TODO: Shared statistics class between reader and writer.
public class BamToFastq implements Closeable
{
    private final BamToFastqConfig mConfig;
    private final List<RegionTaskConsumer> mConsumers;
    private final ReadCacheInterface mReadCache;
    private final PairedFastqWriterInterface mPairedFastqWriter;
    private final Closer mCloser;
    private final ProgressReporter mProgressReporter;

    @Inject
    public BamToFastq(final BamToFastqConfig config, final List<RegionTaskConsumer> consumers, final ReadCacheInterface readCache,
            final PairedFastqWriterInterface pairedFastqWriter, final ProgressReporter progressReporter)
    {
        mConfig = config;
        mConsumers = consumers;
        mReadCache = readCache;
        mCloser = Closer.create();
        mPairedFastqWriter = mCloser.register(pairedFastqWriter);
        mProgressReporter = mCloser.register(progressReporter);
    }

    public void run()
    {
        final long startTimeNanos = System.nanoTime();
        if(mConfig.SilentValidation)
        {
            SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
        }

        mProgressReporter.start();
        if(mConsumers.size() == 1)
        {
            mConsumers.get(0).run();
        }
        else
        {
            for(Thread consumer : mConsumers)
            {
                consumer.start();
            }

            for(Thread consumer : mConsumers)
            {
                try
                {
                    consumer.join();
                }
                catch(InterruptedException e)
                {
                    throw new RuntimeException(e);
                }
            }
        }

        mReadCache.flush();
        if(!mReadCache.isEmpty())
        {
            throw new RuntimeException("There are orphaned reads");
        }
        mReadCache.logStats();
        mPairedFastqWriter.logStats();

        long endTimeNanos = System.nanoTime();
        long durationNanos = endTimeNanos - startTimeNanos;
        float elapsedSeconds = 1.0f * durationNanos / 1_000_000_000;
        BFQ_LOGGER.info(format("BamToFastq complete: elapsedTime(%.3fs) cpuTime(%.3fs)", elapsedSeconds,
                max(1, mConfig.Threads) * elapsedSeconds));
    }

    @Override
    public void close() throws IOException
    {
        mCloser.close();
    }

    public static synchronized void uncaughtExceptionHandler(@Nullable final Thread t, final Throwable e)
    {
        String message = e.getMessage();
        if(t != null)
        {
            message = format("inThread(%s) %s", t.getName(), message);
        }

        Throwable cause = e.getCause();
        if(cause == null)
        {
            BFQ_LOGGER.error(message);
        }
        else
        {
            BFQ_LOGGER.error("{}: {}", message, cause.toString());
        }

        StringWriter stringWriter = new StringWriter();
        e.printStackTrace(new PrintWriter(stringWriter));
        BFQ_LOGGER.debug(stringWriter.toString());
        System.exit(1);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);
        BamToFastqConfig config = new BamToFastqConfig(configBuilder);

        Injector injector = Guice.createInjector(Stage.PRODUCTION, new DefaultModule(config));
        BamToFastq app = injector.getInstance(BamToFastq.class);
        try
        {
            app.run();
        }
        catch(RuntimeException e)
        {
            uncaughtExceptionHandler(null, e);
        }
        finally
        {
            try
            {
                app.close();
            }
            catch(IOException e)
            {
            }
        }
    }
}
