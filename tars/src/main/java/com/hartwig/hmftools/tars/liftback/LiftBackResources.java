package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.tars.liftback.features.GenomicAlignmentScorer;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryConfig;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryMerger;

// A factory rather than a plain holder because htsjdk's indexed FASTA reader is not thread-safe: each worker gets its
// own ref genome handle and its own copy of everything holding one.
public final class LiftBackResources
{
    private final PlacementSelector mPlacementSelector;
    private final EnsemblAnnotationIndex mEnsemblAnnotationIndex;
    private final String mRefGenomeFile;
    private final SupplementaryConfig mSupplementary;
    private final ExcludedRegions mExcludedRegions; // nullable: unmaps primaries that land in these regions post-lift

    public LiftBackResources(
            final PlacementSelector placementSelector, final EnsemblAnnotationIndex annotationIndex, final String refGenomeFile,
            final SupplementaryConfig supplementary, final ExcludedRegions excludedRegions)
    {
        mPlacementSelector = placementSelector;
        mEnsemblAnnotationIndex = annotationIndex;
        mRefGenomeFile = refGenomeFile;
        mSupplementary = supplementary;
        mExcludedRegions = excludedRegions;
    }

    public LiftBackGroupProcessor createProcessor()
    {
        RefGenomeInterface refGenome = openRefGenome();

        return new LiftBackGroupProcessor(
                mPlacementSelector,
                new SupplementaryMerger(mEnsemblAnnotationIndex, refGenome, mSupplementary),
                new OverhangGate(refGenome),
                new GenomicAlignmentScorer(refGenome),
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
