package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;

import com.hartwig.hmftools.common.bam.SamRecordUtils;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;

public record SagaAlignment(
        BwaMemAlignment rawAlignment,
        Cigar cigar,
        int queryLength,
        SagaAssembly sagaAssembly
)
{
    public boolean isForward()
    {
        return !SamRecordUtils.isFlagSet(rawAlignment.getSamFlag(), SAMFlag.READ_REVERSE_STRAND);
    }

    public int queryStart()
    {
        return leftClipLength(cigar);
    }

    public int queryEnd()
    {
        return queryLength - rightClipLength(cigar);
    }

    public int queryAlignLength()
    {
        return queryEnd() - queryStart();
    }

    public int sagaStart()
    {
        return rawAlignment.getRefStart();
    }

    public int sagaEnd()
    {
        return rawAlignment.getRefEnd();
    }

    public int sagaAlignLength()
    {
        return sagaEnd() - sagaStart();
    }

    public int sagaLength()
    {
        return sagaAssembly.assemblyLength();
    }

    public int alignScore()
    {
        return rawAlignment.getAlignerScore();
    }

    // How many more bases could've been aligned on the left?
    public int leftUnaligned()
    {
        return min(queryStart(), sagaStart());
    }

    // How many more bases could've been aligned on the right?
    public int rightUnaligned()
    {
        return min(queryLength - queryEnd(), sagaLength() - sagaEnd());
    }

    public SagaVariant sagaVariant()
    {
        return sagaAssembly.variant();
    }

    @NotNull
    @Override
    public String toString()
    {
        return String.format(
                "SagaAlignment(variant=\"%s\", queryAlign=[%d, %d), sagaAlign=[%d, %d), strand=%s, cigar=%s, alignScore=%d)",
                sagaAssembly().variant(), queryStart(), queryEnd(), sagaStart(), sagaEnd(), isForward() ? "+" : "-", cigar(), alignScore());
    }
}
