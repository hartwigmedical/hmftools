package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.utils.Arrays;

import htsjdk.samtools.CigarElement;

public class VariantReadContext
{
    public final int AlignmentStart;
    public final int AlignmentEnd;

    public final byte[] RefBases;

    // read bases and info
    public final byte[] ReadBases;
    public final List<CigarElement> ReadCigar;
    public final int CoreIndexStart;
    public final int VarReadIndex;
    public final int CoreIndexEnd;

    public final Microhomology Homology;
    public final RepeatInfo MaxRepeat;

    private final SimpleVariant mVariant;
    private final String mReadCigarStr;

    public  VariantReadContext(
            final SimpleVariant variant, final int alignmentStart, final int alignmentEnd, final byte[] refBases,
            final byte[] readBases, final List<CigarElement> readCigar,
            final int coreIndexStart, final int varReadIndex, final int coreIndexEnd,
            final Microhomology homology, final RepeatInfo maxRepeat)
    {
        mVariant = variant;
        AlignmentStart = alignmentStart;
        AlignmentEnd = alignmentEnd;
        RefBases = refBases;
        ReadBases = readBases;
        ReadCigar = readCigar;
        CoreIndexStart = coreIndexStart;
        VarReadIndex = varReadIndex;
        CoreIndexEnd = coreIndexEnd;
        Homology = homology;
        MaxRepeat = maxRepeat;

        mReadCigarStr = CigarUtils.cigarStringFromElements(readCigar);
    }

    // read context methods
    public String ref() { return mVariant.ref(); }
    public String alt() { return mVariant.alt(); }

    public int coreLength() { return CoreIndexEnd - CoreIndexStart + 1; }
    public int leftFlankLength() { return CoreIndexStart; }
    public int rightFlankLength() { return ReadBases.length - CoreIndexEnd - 1; }

    public int refIndex() { return mVariant.Position - AlignmentStart; }

    public String coreStr() { return new String(ReadBases, CoreIndexStart, coreLength()); }
    public String leftFlankStr() { return new String(ReadBases, 0, leftFlankLength()); }
    public String rightFlankStr() { return new String(ReadBases, CoreIndexEnd + 1, rightFlankLength()); }

    public String readBases() { return new String(ReadBases); }
    public String refBases() { return new String(RefBases); }

    public final byte[] trinucleotide()
    {
        int refIndex = refIndex();
        return Arrays.subsetArray(RefBases, refIndex - 1, refIndex + 1);
    }

    public final String trinucleotideStr() { return new String(trinucleotide()); }
    public final String readCigar() { return mReadCigarStr; }

    public String toString()
    {
        return format("%s-%s-%s %s pos(%d-%d) index(%d-%d-%d) repeat(%s) homology(%s)",
                leftFlankStr(), coreStr(), rightFlankStr(), mReadCigarStr, AlignmentStart, AlignmentEnd,
                CoreIndexStart, VarReadIndex, CoreIndexEnd, MaxRepeat != null ? MaxRepeat : "", Homology != null ? Homology : "");
    }
}
