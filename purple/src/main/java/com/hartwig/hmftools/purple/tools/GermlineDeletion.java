package com.hartwig.hmftools.purple.tools;

import com.hartwig.hmftools.common.purple.GermlineStatus;

public class GermlineDeletion
{
    public final String SampleId;
    public final String Chromosome;
    public final int RegionStart;
    public final int RegionEnd;
    public final GermlineStatus Status;
    public final int RegionBafCount;
    public final int RegionDWC;
    public final double RegionObsTumorRatio;
    public final double RegionObsNormalRatio;
    public final double RegionRefNormalisedCN;
    public final double PurpleCN;
    public final int PurpleDWC;
    public final double GcContent;
    public final double MajorAlleleCN;
    public final double MinorAlleleCN;
    public final String GeneName;

    public GermlineDeletion(
            final String sampleId, final String chromosome, final int regionStart, final int regionEnd,
            final GermlineStatus status, final int regionBafCount, final int regionDWC, final double regionObsTumorRatio,
            final double regionObsNormalRatio, final double regionRefNormalisedCN, final double purpleCN, final int purpleDWC,
            final double gcContent, final double majorAlleleCN, final double minorAlleleCN, final String geneName)
    {
        SampleId = sampleId;
        Chromosome = chromosome;
        RegionStart = regionStart;
        RegionEnd = regionEnd;
        Status = status;
        RegionBafCount = regionBafCount;
        RegionDWC = regionDWC;
        RegionObsTumorRatio = regionObsTumorRatio;
        RegionObsNormalRatio = regionObsNormalRatio;
        RegionRefNormalisedCN = regionRefNormalisedCN;
        PurpleCN = purpleCN;
        PurpleDWC = purpleDWC;
        GcContent = gcContent;
        MajorAlleleCN = majorAlleleCN;
        MinorAlleleCN = minorAlleleCN;
        GeneName = geneName;
    }
}
