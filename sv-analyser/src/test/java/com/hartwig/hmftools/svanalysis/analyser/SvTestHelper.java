package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;

import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

public class SvTestHelper
{
    public static SvVarData createTestSv(final String varId, final String chrStart, final String chrEnd,
            long posStart, long posEnd, int orientStart, int orientEnd, StructuralVariantType type)
    {
        return createTestSv(varId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, type,
                1, 1, 1, 1, 1);
    }

    // for convenience
    public static SvVarData createDel(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, 1, -1, DEL,
                1, 1, 1, 1, 1);
    }

    public static SvVarData createDup(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, -1, 1, DUP,
                1, 1, 1, 1, 1);
    }

    public static SvVarData createInv(final String varId, final String chromosome, long posStart, long posEnd, int orientation)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, orientation, orientation, INV,
                1, 1, 1, 1, 1);
    }

    public static SvVarData createSgl(final String varId, final String chromosome, long position, int orientation, boolean isNoneSegment)
    {
        SvVarData var = createTestSv(varId, chromosome, "0", position, -1, orientation, -1, SGL,
                1, 0, 1, 0, 1);

        var.setNoneSegment(isNoneSegment);

        return var;
    }

    public static SvVarData createTestSv(final String varId, final String chrStart, final String chrEnd,
            long posStart, long posEnd, int orientStart, int orientEnd, StructuralVariantType type,
            double cnStart, double cnEnd, double cnChgStart, double cnChgEnd, double ploidy)
    {
        StructuralVariantData svData =
                ImmutableStructuralVariantData.builder()
                        .id(varId)
                        .startChromosome(chrStart)
                        .endChromosome(chrEnd)
                        .startPosition(posStart)
                        .endPosition(posEnd)
                        .startOrientation((byte)orientStart)
                        .endOrientation((byte)orientEnd)
                        .startAF(1.0)
                        .adjustedStartAF(1.0)
                        .adjustedStartCopyNumber(cnStart)
                        .adjustedStartCopyNumberChange(cnChgStart)
                        .endAF(1.0)
                        .adjustedEndAF(1.0)
                        .adjustedEndCopyNumber(cnEnd)
                        .adjustedEndCopyNumberChange(cnChgEnd)
                        .ploidy(ploidy)
                        .type(type)
                        .homology("")
                        .vcfId("")
                        .insertSequence("")
                        .filter("PASS")
                        .imprecise(false)
                        .qualityScore(0.0)
                        .event("")
                        .startTumourVariantFragmentCount(0)
                        .startTumourReferenceFragmentCount(0)
                        .startNormalVariantFragmentCount(0)
                        .startNormalReferenceFragmentCount(0)
                        .endTumourVariantFragmentCount(0)
                        .endTumourReferenceFragmentCount(0)
                        .endNormalVariantFragmentCount(0)
                        .endNormalReferenceFragmentCount(0)
                        .startIntervalOffsetStart(0)
                        .startIntervalOffsetEnd(0)
                        .endIntervalOffsetStart(0)
                        .endIntervalOffsetEnd(0)
                        .inexactHomologyOffsetStart(0)
                        .inexactHomologyOffsetEnd(0)
                        .startLinkedBy("")
                        .endLinkedBy("")
                        .build();

        SvVarData var = new SvVarData(svData);

        String startArm = getChromosomalArm(var.chromosome(true), var.position(true));

        String endArm = "";
        if(!var.isNullBreakend())
            endArm = getChromosomalArm(var.chromosome(false), var.position(false));
        else
            endArm = CHROMOSOME_ARM_P;

        var.setChromosomalArms(startArm, endArm);

        return var;
    }

}
