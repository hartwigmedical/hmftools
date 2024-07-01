package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.quality.ArtefactContext;
import com.hartwig.hmftools.sage.quality.UltimaQualModel;

import htsjdk.samtools.CigarElement;

public class VariantReadContext
{
    public final int AlignmentStart; // alignments of the flanks using read bases
    public final int AlignmentEnd;

    // read bases and info
    public final byte[] ReadBases;
    public final int CoreIndexStart; // in the read to the start of the core
    public final int VarIndex; // index in the read bases of the variant's position
    public final int CoreIndexEnd;

    public final Microhomology Homology;
    public final RepeatInfo MaxRepeat;
    public final List<RepeatInfo> AllRepeats;

    public final int CorePositionStart;
    public final int CorePositionEnd;

    public final byte[] RefBases; // captured for the core

    private final SimpleVariant mVariant;

    private final String mReadCigarStr;

    private ArtefactContext mArtefactContext;
    private UltimaQualModel mUltimaQualModel;

    private RepeatInfo mMaxRefRepeat; // maximum repeat in the reference, only written to the VCF for downstream usage (ie repeat sites)
    private String mExtendedRefBases;

    public VariantReadContext(
            final SimpleVariant variant, final int alignmentStart, final int alignmentEnd, final byte[] refBases,
            final byte[] readBases, final List<CigarElement> readCigar, final int coreIndexStart, final int varIndex, final int coreIndexEnd,
            final Microhomology homology, final RepeatInfo maxRepeat, final List<RepeatInfo> allRepeats,
            final int corePositionStart, final int corePositionEnd)
    {
        mVariant = variant;
        AlignmentStart = alignmentStart;
        AlignmentEnd = alignmentEnd;
        RefBases = refBases;
        ReadBases = readBases;
        CoreIndexStart = coreIndexStart;
        VarIndex = varIndex;
        CoreIndexEnd = coreIndexEnd;
        Homology = homology;
        MaxRepeat = maxRepeat;
        AllRepeats = allRepeats;

        mReadCigarStr = CigarUtils.cigarStringFromElements(readCigar);

        CorePositionStart = corePositionStart;
        CorePositionEnd = corePositionEnd;

        mArtefactContext = null;
        mUltimaQualModel = null;
        mMaxRefRepeat = null;
        mExtendedRefBases = null;
    }

    // read context methods
    public SimpleVariant variant() { return mVariant; }
    public String ref() { return mVariant.ref(); }
    public String alt() { return mVariant.alt(); }

    public boolean hasHomology() { return Homology != null; }
    public int coreLength() { return CoreIndexEnd - CoreIndexStart + 1; }
    public int leftCoreLength() { return VarIndex - CoreIndexStart; }
    public int rightCoreLength() { return CoreIndexEnd - VarIndex; }
    public int leftFlankLength() { return CoreIndexStart; }
    public int rightFlankLength() { return ReadBases.length - CoreIndexEnd - 1; }
    public int leftLength() { return VarIndex; } // distance from position index to first read base
    public int rightLength() { return ReadBases.length - VarIndex; } // distance to last base
    public int totalLength() { return ReadBases.length; }

    public int variantRefIndex() { return mVariant.position() - CorePositionStart; }

    public boolean isValid()
    {
        if(CoreIndexStart <= 0 || CoreIndexEnd >= ReadBases.length - 1)
            return false; // implies no flank

        int minCoreLength = mVariant.isIndel() ? MIN_CORE_DISTANCE * 2 : MIN_CORE_DISTANCE * 2 + 1;
        if(coreLength() < minCoreLength)
            return false;

        if(VarIndex <= CoreIndexStart || CoreIndexEnd <= VarIndex) // invalid var index
            return false;

        if(mVariant.Position <= CorePositionStart || mVariant.Position >= CorePositionEnd)
            return false;

        return true;
    }

    public String coreStr() { return new String(ReadBases, CoreIndexStart, coreLength()); }
    public String leftFlankStr() { return new String(ReadBases, 0, leftFlankLength()); }
    public String rightFlankStr() { return new String(ReadBases, CoreIndexEnd + 1, rightFlankLength()); }

    public String readBases() { return new String(ReadBases); }
    public String refBases() { return new String(RefBases); }

    public String homologyBases() { return Homology != null ? Homology.Bases : ""; }
    public int maxRepeatCount() { return MaxRepeat != null ? MaxRepeat.Count : 0; }

    public final byte[] trinucleotide()
    {
        int refIndex = variantRefIndex();
        return Arrays.subsetArray(RefBases, refIndex - 1, refIndex + 1);
    }

    public final String trinucleotideStr() { return new String(trinucleotide()); }
    public final String readCigar() { return mReadCigarStr; }

    public ArtefactContext artefactContext() { return mArtefactContext; }
    public void setArtefactContext(final ArtefactContext context) { mArtefactContext = context; }

    public UltimaQualModel ultimaQualModel() { return mUltimaQualModel; }
    public void setUltimaQualModel(final UltimaQualModel model) { mUltimaQualModel = model; }

    public RepeatInfo refMaxRepeat() { return mMaxRefRepeat; }
    public void setRefMaxRepeat(final RepeatInfo repeatInfo) { mMaxRefRepeat = repeatInfo; }

    public String extendedRefBases() { return mExtendedRefBases; }
    public void setExtendedRefBases(final String refBases) { mExtendedRefBases = refBases; }

    public String toString()
    {
        return format("%s read(%s-%s-%s %s) pos(%d-%d) index(%d-%d-%d) repeat(%s) homology(%s) ref(%s)",
                mVariant, leftFlankStr(), coreStr(), rightFlankStr(), mReadCigarStr, AlignmentStart, AlignmentEnd,
                CoreIndexStart, VarIndex, CoreIndexEnd, MaxRepeat != null ? MaxRepeat : "", Homology != null ? Homology : "", refBases());
    }

    @VisibleForTesting
    public boolean matches(final VariantReadContext other)
    {
        if((Homology != null) != (other.Homology != null))
            return false;

        if(Homology != null && !Homology.Bases.equals(other.Homology.Bases))
            return false;

        if((MaxRepeat != null) != (other.MaxRepeat != null))
            return false;

        if(MaxRepeat != null && !MaxRepeat.matches(other.MaxRepeat))
            return false;

        return VarIndex == other.VarIndex && AlignmentStart == other.AlignmentStart && AlignmentEnd == other.AlignmentEnd
                && refBases().equals(other.refBases()) && readBases().equals(other.readBases()) && readCigar().equals(other.readCigar())
                && CorePositionStart == other.CorePositionStart && CorePositionEnd == other.CorePositionEnd;
    }
}
