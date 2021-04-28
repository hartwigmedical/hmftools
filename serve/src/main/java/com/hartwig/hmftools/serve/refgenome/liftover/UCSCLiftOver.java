package com.hartwig.hmftools.serve.refgenome.liftover;

import java.io.File;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

public class UCSCLiftOver implements LiftOverAlgo {

    @NotNull
    private final LiftOver liftOver;

    @NotNull
    public static UCSCLiftOver fromChainFile(@NotNull String chainFile) {
        return new UCSCLiftOver(new LiftOver(new File(chainFile)));
    }

    private UCSCLiftOver(@NotNull final LiftOver liftOver) {
        this.liftOver = liftOver;
    }

    @Nullable
    @Override
    public LiftOverResult liftOver(@NotNull final String chromosome, final long position) {
        // UCSC expects hg19 format in case v37 is used.
        String ucscChromosome = RefGenomeFunctions.enforceChromosome(chromosome);
        Interval interval = new Interval(ucscChromosome, (int) position, (int) position);
        Interval lifted = liftOver.liftOver(interval);
        if (lifted == null) {
            return null;
        }

        return ImmutableLiftOverResult.builder().chromosome(lifted.getContig()).position(lifted.getStart()).build();
    }
}
