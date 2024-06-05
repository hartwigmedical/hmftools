package com.hartwig.hmftools.sage.vis;

import java.util.List;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class BaseSeqViewModel
{
    public final Boolean LeftIsForwardStrand;
    public final Boolean RightIsForwardStrand;
    public final int FirstBasePos;
    public final int LastBasePos;

    private final List<BaseViewModel> mBases;
    private final int mPosStart;

    public BaseSeqViewModel(final List<BaseViewModel> bases, int posStart, @Nullable final Boolean leftIsForwardStrand,
            @Nullable final Boolean rightIsForwardStrand)
    {
        mBases = bases;
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
        List<BaseViewModel> bases = Lists.newArrayList();
        for(int i = 0; i < baseStr.length(); ++i)
        {
            bases.add(new BaseViewModel(baseStr.charAt(i)));
        }

        return new BaseSeqViewModel(bases, posStart, null, null);
    }

    public static BaseSeqViewModel fromVariant(final VariantReadContext readContext, final String ref, final String alt)
    {
        String rawBases = readContext.readBases();
        int posStart = readContext.variant().Position - readContext.VarIndex;
        if(ref.length() == alt.length())
        {
            return fromStr(rawBases, posStart);
        }

        // del
        if(ref.length() > alt.length())
        {
            int delLen = ref.length() - alt.length();

            // alt is single char
            List<BaseViewModel> bases = Lists.newArrayList();
            for(int i = 0; i <= readContext.VarIndex; ++i)
            {
                bases.add(new BaseViewModel(rawBases.charAt(i)));
            }

            for(int i = 0; i < delLen; ++i)
            {
                bases.add(BaseViewModel.createDelBase());
            }

            for(int i = readContext.VarIndex + 1; i < readContext.totalLength(); ++i)
            {
                bases.add(new BaseViewModel(rawBases.charAt(i)));
            }

            return new BaseSeqViewModel(bases, posStart, null, null);
        }

        // ins
        int insLen = alt.length() - ref.length();

        // ref is single char
        List<BaseViewModel> bases = Lists.newArrayList();
        for(int i = 0; i <= readContext.VarIndex; ++i)
        {
            bases.add(new BaseViewModel(rawBases.charAt(i)));
        }

        bases.get(bases.size() - 1).incRightInsertCount(insLen);

        for(int i = readContext.VarIndex + insLen + 1; i < readContext.totalLength(); ++i)
        {
            bases.add(new BaseViewModel(rawBases.charAt(i)));
        }

        return new BaseSeqViewModel(bases, posStart, null, null);
    }

    public static BaseSeqViewModel fromRead(final SAMRecord read)
    {
        return fromConsensusFragment(read, null, null);
    }

    public static BaseSeqViewModel fromConsensusFragment(final SAMRecord consensusRead, @Nullable final BaseSeqViewModel first,
            @Nullable final BaseSeqViewModel second)
    {
        int unclippedStart = consensusRead.getUnclippedStart();
        String readString = consensusRead.getReadString();
        byte[] baseQuals = consensusRead.getBaseQualities();

        List<BaseViewModel> indexedBases = Lists.newArrayList();
        int baseIdx = 0;
        for(CigarElement cigarElem : consensusRead.getCigar().getCigarElements())
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
                        indexedBases.add(new BaseViewModel(readString.charAt(baseIdx), baseQuals[baseIdx], false, isOverlapped));
                        baseIdx++;
                        break;
                    case S:
                        indexedBases.add(new BaseViewModel(readString.charAt(baseIdx), baseQuals[baseIdx], true, isOverlapped));
                        baseIdx++;
                        break;
                    case I:
                        if(!indexedBases.isEmpty())
                        {
                            indexedBases.get(indexedBases.size() - 1).incRightInsertCount();
                        }

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
        {
            boolean readNegativeStrandFlag = consensusRead.getReadNegativeStrandFlag();
            return new BaseSeqViewModel(indexedBases, unclippedStart, !readNegativeStrandFlag, !readNegativeStrandFlag);
        }

        boolean leftIsForwardStrand = first.FirstBasePos <= second.FirstBasePos ? first.LeftIsForwardStrand : second.LeftIsForwardStrand;
        boolean rightIsForwardStrand = first.LastBasePos >= second.LastBasePos ? first.RightIsForwardStrand : second.RightIsForwardStrand;
        return new BaseSeqViewModel(indexedBases, unclippedStart, leftIsForwardStrand, rightIsForwardStrand);
    }

    public BaseViewModel getBase(int pos)
    {
        int baseIdx = pos - mPosStart;
        if(baseIdx < 0 || baseIdx >= mBases.size())
        {
            return BaseViewModel.createMissingBase();
        }

        return mBases.get(baseIdx);
    }

    public boolean hasOrientation()
    {
        return LeftIsForwardStrand != null && RightIsForwardStrand != null;
    }
}
