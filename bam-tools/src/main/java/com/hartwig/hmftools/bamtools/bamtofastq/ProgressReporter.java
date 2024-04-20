package com.hartwig.hmftools.bamtools.bamtofastq;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.bamtofastq.BamToFastqConfig.BFQ_LOGGER;

import java.io.Closeable;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.inject.Inject;
import com.hartwig.hmftools.bamtools.bamtofastq.BindingAnnotations.TotalRegionCount;

// TODO NEXT: TEST
public class ProgressReporter implements Closeable
{
    private static final int PERIOD_SECONDS = 5;

    private final int mTotalRegionCount;
    private final AtomicInteger mTotalRegionsProcessedCount;
    private final ScheduledExecutorService mScheduler;
    private final Runnable mReporter = new Runnable()
    {
        @Override
        public void run()
        {
            int totalRegionsProcessedCount = mTotalRegionsProcessedCount.get();
            float progress = 1.0f * totalRegionsProcessedCount / mTotalRegionCount;
            BFQ_LOGGER.info("Processed {} / {} regions progress {}", totalRegionsProcessedCount, mTotalRegionCount, format("%.2f%%",
                    100 * progress));
        }
    };

    private ScheduledFuture<?> mReporterHandle;

    @Inject
    public ProgressReporter(@TotalRegionCount int totalRegionCount, final AtomicInteger totalRegionsProcessedCount)
    {
        mTotalRegionCount = totalRegionCount;
        mTotalRegionsProcessedCount = totalRegionsProcessedCount;
        mScheduler = Executors.newScheduledThreadPool(1);

        mReporterHandle = null;
    }

    public void start()
    {
        mReporterHandle = mScheduler.scheduleAtFixedRate(mReporter, PERIOD_SECONDS, PERIOD_SECONDS, TimeUnit.SECONDS);
    }

    @Override
    public void close()
    {
        mReporterHandle.cancel(false);

        mScheduler.shutdown();
        try
        {
            if(!mScheduler.awaitTermination(60, TimeUnit.SECONDS))
            {
                mScheduler.shutdownNow();
                if(!mScheduler.awaitTermination(60, TimeUnit.SECONDS))
                {
                    BFQ_LOGGER.error("Pool did not terminate");
                }
            }
        }
        catch(InterruptedException ie)
        {
            mScheduler.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }
}
