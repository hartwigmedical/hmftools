package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public interface CombinedRegion extends GenomeRegion {

    @NotNull
    @Override
    default String chromosome() {
        return region().chromosome();
    }

    @Override
    default int start() {
        return region().start();
    }

    @Override
    default int end() {
        return region().end();
    }

    default double tumorCopyNumber() {
        return region().tumorCopyNumber();
    }

    default double tumorBAF() {
        return region().tumorBAF();
    }

    default int bafCount() {
        return region().bafCount();
    }

    boolean isInferredBAF();

    @NotNull
    List<FittedRegion> regions();

    @NotNull
    CopyNumberMethod copyNumberMethod();

    boolean isProcessed();

    void setCopyNumberMethod(@NotNull final CopyNumberMethod copyNumberMethod);

    @NotNull
    FittedRegion region();

    void extend(@NotNull final FittedRegion region);

    void extendWithUnweightedAverage(@NotNull final FittedRegion region);

    void extendWithWeightedAverage(@NotNull final FittedRegion region);

    void setTumorCopyNumber(@NotNull final CopyNumberMethod method, double copyNumber);

    void setInferredTumorBAF(double baf);

    void resetDepthWindowCount();

    @NotNull
    SegmentSupport germlineEndSupport();

    @NotNull
    default SegmentSupport support() {
        return region().support();
    }
}
