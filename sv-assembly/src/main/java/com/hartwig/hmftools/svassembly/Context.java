package com.hartwig.hmftools.svassembly;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.svassembly.assembly.Aligner;
import com.hartwig.hmftools.svassembly.assembly.SupportChecker;
import com.hartwig.hmftools.svassembly.processor.Problem;
import com.hartwig.hmftools.svassembly.sam.CachingSAMSource;
import com.hartwig.hmftools.svassembly.sam.CompositeSAMSource;
import com.hartwig.hmftools.svassembly.sam.DirectSAMSource;
import com.hartwig.hmftools.svassembly.sam.NormalisingSource;
import com.hartwig.hmftools.svassembly.sam.RecordNormaliser;
import com.hartwig.hmftools.svassembly.sam.SAMSource;

import org.jetbrains.annotations.Nullable;

public class Context implements AutoCloseable
{
    public final ExecutorService Executor;
    public final com.hartwig.hmftools.svassembly.assembly.Aligner Aligner;
    @Nullable
    public final com.hartwig.hmftools.svassembly.assembly.Aligner AlternateAligner;
    public final SAMSource SAMSource;
    public final RefGenomeSource ReferenceGenome;
    public final SupportChecker SupportChecker;

    public final SVAConfig Config;

    public final List<Problem> Problems = new ArrayList<>();

    public Context(final ExecutorService executor, final com.hartwig.hmftools.svassembly.assembly.Aligner aligner, @Nullable final com.hartwig.hmftools.svassembly.assembly.Aligner alternateAligner, final SAMSource samSource,
            final RefGenomeSource referenceGenome, final SupportChecker supportChecker, final SVAConfig config)
    {
        Executor = executor;
        Aligner = aligner;
        AlternateAligner = alternateAligner;
        SAMSource = samSource;
        ReferenceGenome = referenceGenome;
        SupportChecker = supportChecker;
        Config = config;
    }

    @Override
    public void close()
    {
        Aligner.close();
        if(AlternateAligner != null)
            AlternateAligner.close();
        SAMSource.close();
    }

    private static String osExtension()
    {
        final String osName = System.getProperty("os.name");
        if (osName.contains("Mac"))
            return ".dylib";
        else if (osName.contains("Win"))
            return ".dll";
        else
            return ".so";
    }

    public static Context create(final SVAConfig config)
    {
        final var props = System.getProperties();
        final String candidateBWAPath = "libbwa." + props.getProperty("os.arch") + osExtension();
        if (System.getProperty("LIBBWA_PATH") == null && new File(candidateBWAPath).exists())
            System.setProperty("LIBBWA_PATH", new File(candidateBWAPath).getAbsolutePath());

        final com.hartwig.hmftools.svassembly.assembly.Aligner aligner = new Aligner(config, config.referenceGenomeIndex());
        final com.hartwig.hmftools.svassembly.assembly.Aligner alternateAligner = config.altReferenceGenomeIndex() == null
                ? null
                : new Aligner(config, Objects.requireNonNull(config.altReferenceGenomeIndex()));

        final RefGenomeSource refGenomeSource = RefGenomeSource.loadRefGenome(config.referenceGenomeFile().getAbsolutePath());
        final SupportChecker supportChecker = new SupportChecker(config);

        final RecordNormaliser normaliser = new RecordNormaliser(refGenomeSource, config);

        final List<SAMSource> inputFiles = new ArrayList<>();
        inputFiles.add(openFile(config.referenceGenomeFile(), normaliser, config.bamFile(), "tumor"));
        if(config.germlineBAMFile() != null)
            inputFiles.add(openFile(config.referenceGenomeFile(), normaliser, config.germlineBAMFile(), "germline"));

        final SAMSource innerSAMSource = inputFiles.size() == 1 ? inputFiles.get(0) : new CompositeSAMSource(inputFiles);
        final SAMSource samSource = new CachingSAMSource(innerSAMSource);

        final ExecutorService executor = config.threads() == -1
                ? ForkJoinPool.commonPool()
                : Executors.newFixedThreadPool(Math.max(1, config.threads()), runnable ->
                {
                    final Thread thread = new Thread(runnable);
                    thread.setDaemon(true);
                    return thread;
                });
        return new Context(executor, aligner, alternateAligner, samSource, refGenomeSource, supportChecker, config);
    }

    private static SAMSource openFile(final File referenceGenome, final RecordNormaliser normaliser, final File bamFile, final String tag)
    {
        final DirectSAMSource directSource = new DirectSAMSource(bamFile, referenceGenome, tag);
        return new NormalisingSource(directSource, normaliser);
    }
}
