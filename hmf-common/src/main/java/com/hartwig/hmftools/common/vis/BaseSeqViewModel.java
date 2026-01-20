package com.hartwig.hmftools.common.vis;

import static com.hartwig.hmftools.common.vis.BaseViewModel.MISSING_BASEQ;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class BaseSeqViewModel
{
    public final String Chromosome;
    public final Boolean LeftIsForwardStrand;
    public final Boolean RightIsForwardStrand;
    public final int FirstBasePos;
    public final int LastBasePos;

    private final List<BaseViewModel> mBases;
    private final int mPosStart;

    public BaseSeqViewModel(final List<BaseViewModel> bases, int posStart, @Nullable final Boolean leftIsForwardStrand,
            @Nullable final Boolean rightIsForwardStrand)
    {
        this(bases, null, posStart, leftIsForwardStrand, rightIsForwardStrand);
    }

    private BaseSeqViewModel(final List<BaseViewModel> bases, @Nullable final String chromosome, int posStart,
            @Nullable final Boolean leftIsForwardStrand, @Nullable final Boolean rightIsForwardStrand)
    {
        mBases = bases;
        Chromosome = chromosome;
        mPosStart = posStart;
        LeftIsForwardStrand = leftIsForwardStrand;
        RightIsForwardStrand = rightIsForwardStrand;

        int firstBasePos = -1;
        int lastBasePos = -1;
        for(int i = 0; i < bases.size(); i++)
        {
            if(bases.get(i).isMissing())
            {
                continue;
            }

            if(firstBasePos == -1)
            {
                firstBasePos = i + posStart;
            }

            lastBasePos = i + posStart;
        }

        FirstBasePos = firstBasePos;
        LastBasePos = lastBasePos;
    }

    public static BaseSeqViewModel fromStr(final String baseStr, int posStart)
    {
        return fromStr(baseStr, null, posStart);
    }

    public static BaseSeqViewModel fromStr(final String baseStr, @Nullable final String chromosome, int posStart)
    {
        List<BaseViewModel> bases = Lists.newArrayList();
        for(int i = 0; i < baseStr.length(); ++i)
        {
            bases.add(new BaseViewModel(baseStr.charAt(i)));
        }

        return new BaseSeqViewModel(bases, chromosome, posStart, null, null);
    }

    public static BaseSeqViewModel fromRead(final SAMRecord read)
    {
        return fromRead(read, null);
    }

    public static BaseSeqViewModel fromRead(final SAMRecord read, @Nullable final Integer unclippedStartOverride)
    {
        return fromConsensusFragment(read, null, null, unclippedStartOverride);
    }

    public static BaseSeqViewModel fromConsensusFragment(final SAMRecord consensusRead, @Nullable final BaseSeqViewModel first,
            @Nullable final BaseSeqViewModel second)
    {
        return fromConsensusFragment(consensusRead, first, second, null);
    }

    private static BaseSeqViewModel fromConsensusFragment(final SAMRecord consensusRead, @Nullable final BaseSeqViewModel first,
            @Nullable final BaseSeqViewModel second, @Nullable final Integer unclippedStartOverride)
    {
        int unclippedStart = unclippedStartOverride == null ? consensusRead.getUnclippedStart() : unclippedStartOverride;
        List<CigarElement> cigarElements = consensusRead.getCigar().getCigarElements();
        byte[] bases = consensusRead.getReadBases();
        byte[] baseQuals = consensusRead.getBaseQualities();
        boolean readNegativeStrandFlag = consensusRead.getReadNegativeStrandFlag();
        return create(unclippedStart, cigarElements, bases, baseQuals, readNegativeStrandFlag, first, second);
    }

    public static BaseSeqViewModel create(int unclippedStart, final List<CigarElement> cigarElements, final byte[] bases,
            final byte[] baseQuals, boolean readNegativeStrandFlag)
    {
        return create(unclippedStart, cigarElements, bases, baseQuals, readNegativeStrandFlag, null, null);
    }

    public static BaseSeqViewModel create(int unclippedStart, final List<CigarElement> cigarElements, final byte[] bases,
            final byte[] baseQuals, boolean readNegativeStrandFlag, @Nullable final BaseSeqViewModel first,
            @Nullable final BaseSeqViewModel second)
    {
        List<BaseViewModel> indexedBases = Lists.newArrayList();
        int baseIdx = 0;
        for(CigarElement cigarElem : cigarElements)
        {
            CigarOperator cigarOp = cigarElem.getOperator();
            int elemLen = cigarElem.getLength();
            for(int i = 0; i < elemLen; i++)
            {
                BaseViewModel firstBase = first == null ? BaseViewModel.createMissingBase() : first.getBase(baseIdx + unclippedStart);
                BaseViewModel secondBase = second == null ? BaseViewModel.createMissingBase() : second.getBase(baseIdx + unclippedStart);
                boolean isOverlapped = !firstBase.isMissing() && !secondBase.isMissing();

                switch(cigarOp)
                {
                    case M:
                    case EQ:
                    case X:
                        indexedBases.add(new BaseViewModel((char) bases[baseIdx], baseQuals[baseIdx], false, isOverlapped));
                        baseIdx++;
                        break;
                    case S:
                        indexedBases.add(new BaseViewModel((char) bases[baseIdx], baseQuals[baseIdx], true, isOverlapped));
                        baseIdx++;
                        break;
                    case I:
                        if(!indexedBases.isEmpty())
                            indexedBases.get(indexedBases.size() - 1).incRightInsertCount((char) bases[baseIdx], baseQuals[baseIdx]);

                        baseIdx++;
                        break;
                    case D:
                        indexedBases.add(BaseViewModel.createDelBase(isOverlapped));
                        break;
                    case H:
                    case N:
                        indexedBases.add(BaseViewModel.createMissingBase());
                        break;
                }
            }
        }

        if(first == null || second == null)
            return new BaseSeqViewModel(indexedBases, unclippedStart, !readNegativeStrandFlag, !readNegativeStrandFlag);

        boolean leftIsForwardStrand = first.FirstBasePos <= second.FirstBasePos ? first.LeftIsForwardStrand : second.LeftIsForwardStrand;
        boolean rightIsForwardStrand = first.LastBasePos >= second.LastBasePos ? first.RightIsForwardStrand : second.RightIsForwardStrand;
        return new BaseSeqViewModel(indexedBases, unclippedStart, leftIsForwardStrand, rightIsForwardStrand);
    }

    public static BaseSeqViewModel fromStringWithCigar(final String baseStr, final List<CigarElement> cigarElements, int posStart)
    {
        List<BaseViewModel> indexedBases = Lists.newArrayList();
        int baseIdx = 0;
        for(CigarElement cigarElem : cigarElements)
        {
            CigarOperator cigarOp = cigarElem.getOperator();
            int elemLen = cigarElem.getLength();
            for(int i = 0; i < elemLen; i++)
            {
                switch(cigarOp)
                {
                    case M:
                    case EQ:
                    case X:
                        indexedBases.add(new BaseViewModel(baseStr.charAt(baseIdx)));
                        baseIdx++;
                        break;
                    case S:
                        indexedBases.add(new BaseViewModel(baseStr.charAt(baseIdx), MISSING_BASEQ, true, false));
                        baseIdx++;
                        break;
                    case I:
                        if(!indexedBases.isEmpty())
                            indexedBases.get(indexedBases.size() - 1).incRightInsertCount(baseStr.charAt(baseIdx), MISSING_BASEQ);

                        baseIdx++;
                        break;
                    case D:
                        indexedBases.add(BaseViewModel.createDelBase(false));
                        break;
                    case H:
                    case N:
                        indexedBases.add(BaseViewModel.createMissingBase());
                        break;
                }
            }
        }

        return new BaseSeqViewModel(indexedBases, posStart, null, null);
    }

    public BaseViewModel getBase(int pos)
    {
        int baseIdx = pos - mPosStart;
        if(baseIdx < 0 || baseIdx >= mBases.size())
            return BaseViewModel.createMissingBase();

        return mBases.get(baseIdx);
    }

    public boolean hasOrientation()
    {
        return LeftIsForwardStrand != null && RightIsForwardStrand != null;
    }

    public int mismatchCount(final BaseSeqViewModel refBases)
    {
        int mismatches = 0;
        for(int i = refBases.FirstBasePos; i <= refBases.LastBasePos; ++i)
        {
            BaseViewModel refBase = refBases.getBase(i);
            BaseViewModel base = getBase(i);

            if(!base.isSoftClip() && base.hasCharBase() && refBase.hasCharBase() && refBase.charBase() != base.charBase())
                mismatches++;
        }

        return mismatches;
    }

    public void clearSoftClips()
    {
        mBases.forEach(BaseViewModel::clearSoftClip);
    }
}
