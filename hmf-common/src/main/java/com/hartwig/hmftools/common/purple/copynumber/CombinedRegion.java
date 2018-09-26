package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public interface CombinedRegion extends GenomeRegion {

    @NotNull
    @Override
    default String chromosome() {
        return region().chromosome();
    }

    @Override
    default long start() {
        return region().start();
    }

    @Override
    default long end() {
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

    List<FittedRegion> regions();

    CopyNumberMethod copyNumberMethod();

    boolean isProcessed();

    void setCopyNumberMethod(CopyNumberMethod copyNumberMethod);

    FittedRegion region();

    void extend(final FittedRegion region);

    void extendWithUnweightedAverage(final FittedRegion region);

    void extendWithBAFWeightedAverage(final FittedRegion region);

    void setTumorCopyNumber(@NotNull final CopyNumberMethod method, double copyNumber);

    void setInferredTumorBAF(double baf);

}
