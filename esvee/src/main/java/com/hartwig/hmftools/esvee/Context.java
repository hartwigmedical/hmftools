package com.hartwig.hmftools.esvee;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.esvee.assembly.Aligner;
import com.hartwig.hmftools.esvee.assembly.SupportChecker;
import com.hartwig.hmftools.esvee.processor.Problem;
import com.hartwig.hmftools.esvee.read.CachingSAMSource;
import com.hartwig.hmftools.esvee.read.CompositeSAMSource;
import com.hartwig.hmftools.esvee.read.DirectSAMSource;
import com.hartwig.hmftools.esvee.read.RecordNormaliser;
import com.hartwig.hmftools.esvee.read.SAMSource;

public class Context implements AutoCloseable
{
    public final ExecutorService Executor;
    public final com.hartwig.hmftools.esvee.assembly.Aligner Aligner;

    public final SAMSource SAMSource;
    public final RefGenomeSource ReferenceGenome;
    public final SupportChecker SupportChecker;

    public final SvConfig Config;

    public final List<Problem> Problems = new ArrayList<>();

    public Context(
            final ExecutorService executor, final Aligner aligner, final SAMSource samSource,
            final RefGenomeSource referenceGenome, final SupportChecker supportChecker, final SvConfig config)
    {
        Executor = executor;
        Aligner = aligner;
        SAMSource = samSource;
        ReferenceGenome = referenceGenome;
        SupportChecker = supportChecker;
        Config = config;
    }

    @Override
    public void close()
    {
        Aligner.close();
        SAMSource.close();
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
        final SupportChecker supportChecker = new SupportChecker();

        final RecordNormaliser normaliser = new RecordNormaliser(refGenomeSource);

        List<SAMSource> inputBams = Lists.newArrayList();

        for(int i = 0; i < config.SampleNames.size(); ++i)
        {
            String sampleId = config.SampleNames.get(i);
            String bamFile = config.BamFiles.get(i);

            inputBams.add(openFile(config.RefGenomeFile, normaliser, bamFile, sampleId));
        }

        final SAMSource innerSAMSource = inputBams.size() == 1 ? inputBams.get(0) : new CompositeSAMSource(inputBams);
        final SAMSource samSource = new CachingSAMSource(innerSAMSource);

        final ExecutorService executor = config.Threads == -1
                ? ForkJoinPool.commonPool()
                : Executors.newFixedThreadPool(Math.max(1, config.Threads), runnable ->
                {
                    final Thread thread = new Thread(runnable);
                    thread.setDaemon(true);
                    return thread;
                });

        return new Context(executor, aligner, samSource, refGenomeSource, supportChecker, config);
    }

    private static SAMSource openFile(
            final String referenceGenome, final RecordNormaliser normaliser, final String bamFile, final String sampleId)
    {
        return new DirectSAMSource(new File(bamFile), new File(referenceGenome), sampleId);
    }
}
