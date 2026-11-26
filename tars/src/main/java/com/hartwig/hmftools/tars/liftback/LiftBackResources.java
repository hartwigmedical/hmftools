package com.hartwig.hmftools.tars.liftback;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.tars.liftback.merge.SupplementaryConfig;
import com.hartwig.hmftools.tars.liftback.merge.SupplementaryMerger;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;
import com.hartwig.hmftools.tars.liftback.select.AlignmentSelector;

// A factory rather than a plain holder because htsjdk's indexed FASTA reader is not thread-safe: each worker gets its
// own ref genome handle and its own copy of everything holding one.
public final class LiftBackResources
{
    private final AlignmentSelector mAlignmentSelector;
    private final EnsemblAnnotationIndex mEnsemblAnnotationIndex;
    private final String mRefGenomeFile;
    private final SupplementaryConfig mSupplementary;
    private final ExcludedRegions mExcludedRegions; // nullable: unmaps primaries that land in these regions post-lift

    public LiftBackResources(
            final AlignmentSelector alignmentSelector, final EnsemblAnnotationIndex annotationIndex,
            final String refGenomeFile,
            final SupplementaryConfig supplementary, final ExcludedRegions excludedRegions)
    {
        mAlignmentSelector = alignmentSelector;
        mEnsemblAnnotationIndex = annotationIndex;
        mRefGenomeFile = refGenomeFile;
        mSupplementary = supplementary;
        mExcludedRegions = excludedRegions;
    }

    public LiftBackGroupProcessor createProcessor()
    {
        RefGenomeSource refGenome = RefGenomeSource.loadRefGenome(mRefGenomeFile);
        if(refGenome == null)
        {
            throw new IllegalStateException(String.format("failed to load ref genome: %s", mRefGenomeFile));
        }

        return new LiftBackGroupProcessor(
                mAlignmentSelector,
                new SupplementaryMerger(mEnsemblAnnotationIndex, refGenome, mSupplementary),
                new OverhangGate(refGenome),
                new GenomicAlignmentScorer(refGenome),
                refGenome,
                mExcludedRegions);
    }
}
