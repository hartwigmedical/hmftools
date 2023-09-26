package com.hartwig.hmftools.sieve.annotate;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class Util
{
    // TODO(m_cooper): Use the original implementation.
    public static boolean isNotProperReadPair(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
        {
            return false;
        }

        // or a fragment length outside the observed distribution
        // TODO(m_cooper): Make this configurable.
        if(abs(read.getInferredInsertSize()) > 1000)
        {
            return true;
        }

        // an unmapped mate
        if(mateUnmapped(read))
        {
            return true;
        }

        if(read.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!read.getReferenceName().equals(read.getMateReferenceName()))
            {
                return true;
            }

            // inversion
            return read.getReadNegativeStrandFlag() == mateNegativeStrand(read);
        }

        return false;
    }

    public static boolean isSoftClipped(final SAMRecord read)
    {
        for(CigarElement el : read.getCigar().getCigarElements())
        {
            final CigarOperator op = el.getOperator();
            if(op == CigarOperator.SOFT_CLIP)
            {
                return true;
            }
        }

        return false;
    }

    @Nullable
    public static Integer getNM(final SAMRecord read)
    {
        Object rawAttr = read.getAttribute("NM");
        return rawAttr == null ? null : (int) rawAttr;
    }
}
