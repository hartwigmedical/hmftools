package com.hartwig.hmftools.datamodel.purple;

public final class GermlineDeletion {
    public final String GeneName;
    public final String Chromosome;
    public final String ChromosomeBand;
    public final int RegionStart;
    public final int RegionEnd;
    public final int DepthWindowCount;
    public final int ExonStart;
    public final int ExonEnd;
    public final GermlineDetectionMethod DetectionMethod;
    public final GermlineStatus NormalStatus;
    public final GermlineStatus TumorStatus;
    public final double GermlineCopyNumber;
    public final double TumorCopyNumber;
    public final String Filter;
    public final int CohortFrequency;
    public final boolean Reported;

    public GermlineDeletion(
            final String geneName, final String chromosome, final String chromosomeBand, final int regionStart, final int regionEnd,
            final int depthWindowCount, final int exonStart, final int exonEnd, final GermlineDetectionMethod detectionMethod,
            final GermlineStatus normalStatus, final GermlineStatus tumorStatus, final double germlineCopyNumber, final double tumorCopyNumber,
            final String filter, final int cohortFrequency, final boolean reported) {
        GeneName = geneName;
        Chromosome = chromosome;
        ChromosomeBand = chromosomeBand;
        RegionStart = regionStart;
        RegionEnd = regionEnd;
        DepthWindowCount = depthWindowCount;
        ExonStart = exonStart;
        ExonEnd = exonEnd;
        DetectionMethod = detectionMethod;
        NormalStatus = normalStatus;
        TumorStatus = tumorStatus;
        TumorCopyNumber = tumorCopyNumber;
        GermlineCopyNumber = germlineCopyNumber;
        Filter = filter;
        CohortFrequency = cohortFrequency;
        Reported = reported;
    }
}
