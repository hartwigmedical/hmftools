package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;
import com.hartwig.hmftools.tars.liftback.supplementary.AnnotatedJunctionIndex;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryConfig;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryResolver;

// The read-only inputs every LiftBackWorker shares, plus the per-worker assembly that turns them into a processor.
// The ref genome is why this is a factory rather than a plain holder: htsjdk's indexed FASTA reader is not
// thread-safe, so each worker gets its own handle and its own copy of everything that keeps a reference to one.
public final class LiftBackResources
{
    private final LiftBackDiscriminator mDiscriminator;
    private final AnnotatedJunctionIndex mJunctionIndex;
    private final String mRefGenomeFile;
    private final SupplementaryConfig mSupplementary;
    private final ExcludedRegions mExcludedRegions; // nullable: drop fragments here before lifting

    public LiftBackResources(
            final LiftBackDiscriminator discriminator, final AnnotatedJunctionIndex junctionIndex, final String refGenomeFile,
            final SupplementaryConfig supplementary, final ExcludedRegions excludedRegions)
    {
        mDiscriminator = discriminator;
        mJunctionIndex = junctionIndex;
        mRefGenomeFile = refGenomeFile;
        mSupplementary = supplementary;
        mExcludedRegions = excludedRegions;
    }

    // One call per worker. The discriminator, junction index and excluded regions are immutable and shared; the ref
    // genome handle and the two engines that hold on to it are not, so they are built fresh here.
    public LiftBackGroupProcessor createProcessor()
    {
        RefGenomeInterface refGenome = openRefGenome();

        return new LiftBackGroupProcessor(
                mDiscriminator,
                new SupplementaryResolver(mJunctionIndex, refGenome, mSupplementary),
                new OverhangGate(refGenome),
                refGenome,
                mExcludedRegions);
    }

    private RefGenomeInterface openRefGenome()
    {
        RefGenomeSource refGenome = RefGenomeSource.loadRefGenome(mRefGenomeFile);

        if(refGenome == null)
        {
            TARS_LOGGER.error("failed to load ref genome: {}", mRefGenomeFile);
            System.exit(1);
        }

        return refGenome;
    }
}
