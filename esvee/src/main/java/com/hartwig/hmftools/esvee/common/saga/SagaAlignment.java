package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarAlignedLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.cigarBaseLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.getReadIndexFromPosition;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;

import java.util.List;

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
    public void validate()
    {
        if(!(queryStart() >= 0 && queryEnd() <= queryLength))
        {
            throw new IllegalArgumentException();
        }
        if(!(sagaStart() >= 0 && sagaEnd() <= sagaLength()))
        {
            throw new IllegalArgumentException();
        }
        if(!(cigarBaseLength(cigar) == queryLength))
        {
            throw new IllegalArgumentException();
        }
        if(!(cigarAlignedLength(cigar) <= sagaLength()))
        {
            throw new IllegalArgumentException();
        }
    }

    public boolean isForward()
    {
        return !SamRecordUtils.isFlagSet(rawAlignment.getSamFlag(), SAMFlag.READ_REVERSE_STRAND);
    }

    public int queryStart()
    {
        return isForward() ? leftClipLength(cigar) : rightClipLength(cigar);
    }

    public int queryEnd()
    {
        return queryLength - (isForward() ? rightClipLength(cigar) : leftClipLength(cigar));
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
        return sagaAssembly.length();
    }

    public int alignScore()
    {
        return rawAlignment.getAlignerScore();
    }

    // How many more bases could've been aligned at the start of the sequences?
    public int startUnaligned()
    {
        return min(queryStart(), sagaStart());
    }

    // How many more bases could've been aligned at the end of the sequences?
    public int endUnaligned()
    {
        return min(queryLength - queryEnd(), sagaLength() - sagaEnd());
    }

    public SagaVariant sagaVariant()
    {
        return sagaAssembly.variant();
    }

    public List<Integer> queryJunctionOffsets()
    {
        // The indices of the SAGA junctions, mapped into the query sequence range, with extrapolation if required.
        return sagaAssembly().junctionOffsets().stream()
                .map(this::sagaIndexToQueryIndex)
                .toList();
    }

    int sagaIndexToQueryIndex(int sagaIndex)
    {
        // Extrapolate outside of the aligned range.
        if(sagaIndex <= sagaStart())
        {
            return queryStart() - (sagaStart() - sagaIndex);
        }
        if(sagaIndex >= sagaEnd())
        {
            return queryEnd() + (sagaIndex - sagaEnd());
        }
        // On a D element, use the next index.
        int result = getReadIndexFromPosition(sagaStart(), cigar.getCigarElements(), sagaIndex, 1, false);
        if(!(result >= 0 && result <= queryLength))
        {
            // Should be impossible because the out of bounds cases were already handled.
            throw new RuntimeException();
        }
        return result;
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
