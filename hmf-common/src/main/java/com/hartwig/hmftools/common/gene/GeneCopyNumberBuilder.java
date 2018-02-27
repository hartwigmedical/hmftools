package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.hmfslicer.HmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.zipper.RegionZipperHandler;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class GeneCopyNumberBuilder implements RegionZipperHandler<PurpleCopyNumber, HmfExonRegion> {

    private final ImmutableGeneCopyNumber.Builder builder;

    private double minCopyNumber = Double.MAX_VALUE;
    private double minMinorAllelePloidy = Double.MAX_VALUE;
    private double maxCopyNumber = -Double.MAX_VALUE;
    private int somaticCount;
    private int homCount;
    private int het2HomCount;

    @Nullable
    private PurpleCopyNumber previous;

    private double previousCopyNumber = -Double.MAX_VALUE;
    private HmfExonRegion exon;
    private PurpleCopyNumber copyNumber;

    private int minRegions = 0;
    private long minRegionStart = 0;
    private long minRegionEnd = 0;
    private SegmentSupport minRegionStartSupport = SegmentSupport.NONE;
    private SegmentSupport minRegionEndSupport = SegmentSupport.NONE;
    private CopyNumberMethod minRegionMethod = CopyNumberMethod.UNKNOWN;

    private final String gene;
    private int nonsenseBiallelicCount;
    private int nonsenseNonBiallelicCount;
    private double nonsenseNonBiallelicPloidy;
    private int spliceBiallelicCount;
    private int spliceNonBiallelicCount;
    private double spliceNonBiallelicPloidy;
    private int missenseBiallelicCount;
    private int missenseNonBiallelicCount;
    private double missenseNonBiallelicPloidy;

    GeneCopyNumberBuilder(@NotNull final HmfGenomeRegion gene) {
        this.gene = gene.gene();
        builder = ImmutableGeneCopyNumber.builder()
                .from(gene)
                .minRegionStart(gene.start())
                .minRegionEnd(gene.end())
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE);
    }

    @Override
    public void enterChromosome(@NotNull final String chromosome) {
        // IGNORE
    }

    @Override
    public void primary(@NotNull final PurpleCopyNumber copyNumber) {
        this.copyNumber = copyNumber;
        if (exon != null) {
            addOverlap(exon, this.copyNumber);
        }
    }

    @Override
    public void secondary(@NotNull final HmfExonRegion exon) {
        this.exon = exon;
        if (copyNumber != null) {
            addOverlap(this.exon, copyNumber);
        }
    }

    public void somatic(@NotNull final PurityAdjustedSomaticVariant variant) {
        if (variant.gene().equals(gene)) {
            if (variant.worstCodingEffect().equals(CodingEffect.NONSENSE_OR_FRAMESHIFT.toString())) {
                if (variant.biallelic()) {
                    nonsenseBiallelicCount++;
                } else {
                    nonsenseNonBiallelicCount++;
                    nonsenseNonBiallelicPloidy += variant.ploidy();
                }
            } else if (variant.worstCodingEffect().equals(CodingEffect.SPLICE.toString())) {
                if (variant.biallelic()) {
                    spliceBiallelicCount++;
                } else {
                    spliceNonBiallelicCount++;
                    spliceNonBiallelicPloidy += variant.ploidy();
                }
            } else if (variant.worstCodingEffect().equals(CodingEffect.MISSENSE.toString())) {
                if (variant.biallelic()) {
                    missenseBiallelicCount++;
                } else {
                    missenseNonBiallelicCount++;
                    missenseNonBiallelicPloidy += variant.ploidy();
                }
            }

        }
    }

    private void addOverlap(@NotNull final HmfExonRegion exon, @NotNull final PurpleCopyNumber copyNumber) {
        long overlap = exon.overlappingBases(copyNumber);
        if (overlap > 0) {
            double currentCopyNumber = copyNumber.averageTumorCopyNumber();

            maxCopyNumber = Math.max(maxCopyNumber, currentCopyNumber);
            minMinorAllelePloidy = Math.min(minMinorAllelePloidy, copyNumber.minorAllelePloidy());

            if (!Doubles.equal(currentCopyNumber, previousCopyNumber)) {
                switch (copyNumber.method()) {
                    case GERMLINE_HOM_DELETION:
                        homCount++;
                        break;
                    case GERMLINE_HET2HOM_DELETION:
                        het2HomCount++;
                        break;
                    default:
                        somaticCount++;
                }
            }

            if (isUnprocessedCopyNumberRegion(copyNumber)) {
                if (Doubles.lessThan(currentCopyNumber, minCopyNumber)) {
                    minRegions = 1;
                    minCopyNumber = currentCopyNumber;
                    minRegionStart = copyNumber.start();
                    minRegionStartSupport = copyNumber.segmentStartSupport();
                    minRegionEnd = copyNumber.end();
                    minRegionEndSupport = copyNumber.segmentEndSupport();
                    minRegionMethod = copyNumber.method();

                } else if (Doubles.equal(currentCopyNumber, minCopyNumber)) {
                    minRegionEnd = copyNumber.end();
                    minRegionEndSupport = copyNumber.segmentEndSupport();
                    minRegionMethod = copyNumber.method();

                    if (!Doubles.equal(currentCopyNumber, previousCopyNumber)) {
                        minRegions++;
                    }
                }
            }

            previousCopyNumber = currentCopyNumber;
            previous = copyNumber;
        }
    }

    private boolean isUnprocessedCopyNumberRegion(@NotNull final PurpleCopyNumber copyNumber) {
        return previous == null || !previous.equals(copyNumber);
    }

    @NotNull
    public GeneCopyNumber build() {
        return builder.maxCopyNumber(maxCopyNumber)
                .minRegionStartSupport(minRegionStartSupport)
                .minRegionEndSupport(minRegionEndSupport)
                .minRegionMethod(minRegionMethod)
                .minRegionStart(minRegionStart)
                .minRegionEnd(minRegionEnd)
                .minCopyNumber(minCopyNumber)
                .somaticRegions(somaticCount)
                .germlineHomRegions(homCount)
                .germlineHet2HomRegions(het2HomCount)
                .minRegions(minRegions)
                .nonsenseBiallelicCount(nonsenseBiallelicCount)
                .nonsenseNonBiallelicCount(nonsenseNonBiallelicCount)
                .nonsenseNonBiallelicPloidy(nonsenseNonBiallelicPloidy)
                .spliceBiallelicCount(spliceBiallelicCount)
                .spliceNonBiallelicCount(spliceNonBiallelicCount)
                .spliceNonBiallelicPloidy(spliceNonBiallelicPloidy)
                .missenseBiallelicCount(missenseBiallelicCount)
                .missenseNonBiallelicCount(missenseNonBiallelicCount)
                .missenseNonBiallelicPloidy(missenseNonBiallelicPloidy)
                .minMinorAllelePloidy(minMinorAllelePloidy)
                .build();
    }
}
