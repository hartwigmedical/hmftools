package com.hartwig.hmftools.esvee;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.esvee.assembly.Aligner;
import com.hartwig.hmftools.esvee.assembly.SupportChecker;
import com.hartwig.hmftools.esvee.processor.Problem;

public class Context implements AutoCloseable
{
    public final ExecutorService Executor;
    public final com.hartwig.hmftools.esvee.assembly.Aligner Aligner;

    public final RefGenomeSource ReferenceGenome;

    public final SvConfig Config;

    public final List<Problem> Problems = new ArrayList<>();

    public Context(
            final ExecutorService executor, final Aligner aligner,
            final RefGenomeSource referenceGenome, final SvConfig config)
    {
        Executor = executor;
        Aligner = aligner;
        ReferenceGenome = referenceGenome;
        Config = config;
    }

    @Override
    public void close()
    {
        Aligner.close();
    }

    private static String osExtension()
    {
        final String osName = System.getProperty("os.name");
        if(osName.contains("Mac"))
            return ".dylib";
        else if(osName.contains("Win"))
            return ".dll";
        else
            return ".so";
    }

    public static Context create(final SvConfig config)
    {
        final var props = System.getProperties();
        final String candidateBWAPath = "libbwa." + props.getProperty("os.arch") + osExtension();
        if(System.getProperty("LIBBWA_PATH") == null && new File(candidateBWAPath).exists())
            System.setProperty("LIBBWA_PATH", new File(candidateBWAPath).getAbsolutePath());

        Aligner aligner = new Aligner(config, new File(config.RefGenomeImageFile));

        final RefGenomeSource refGenomeSource = RefGenomeSource.loadRefGenome(config.RefGenomeFile);

        final ExecutorService executor = config.Threads == -1
                ? ForkJoinPool.commonPool()
                : Executors.newFixedThreadPool(Math.max(1, config.Threads), runnable ->
                {
                    final Thread thread = new Thread(runnable);
                    thread.setDaemon(true);
                    return thread;
                });

        return new Context(executor, aligner, refGenomeSource, config);
    }
}
