package com.hartwig.hmftools.svanalysis.types;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

public class CopyNumberNoneSegment
{
    public int Id;
    public final String SampleId;
    public final String Chromosome;
    public final long Position;
    public final byte Orientation;
    public final StructuralVariantType Type;
    public final double CopyNumber;
    public final double CopyNumberChange;
    public final double Ploidy;

    public CopyNumberNoneSegment(
            final String sampleId,
            int id,
            final String chromosome,
            final long position,
            final byte orientation,
            final StructuralVariantType type,
            final double copyNumber,
            final double copyNumberChange,
            final double ploidy)
    {
        SampleId = sampleId;
        Id = id;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        Type = StructuralVariantType.SGL;
        CopyNumber = copyNumber;
        CopyNumberChange = copyNumberChange;
        Ploidy = ploidy;
    }

    public CopyNumberNoneSegment(final String sampleId, int id, final PurpleCopyNumber purpleRecord)
    {
        Id = id;
        SampleId = sampleId;
        Chromosome = purpleRecord.chromosome();
        Position = purpleRecord.start();
        Orientation = 1;
        Type = StructuralVariantType.SGL;
        CopyNumber = purpleRecord.averageTumorCopyNumber();
        CopyNumberChange = 0;
        Ploidy = 1;
    }

}
