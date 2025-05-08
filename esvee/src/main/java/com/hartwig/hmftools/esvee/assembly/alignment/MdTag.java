package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;

public class MdTag
{
    private final String mTag;
    private final List<MdTagElement> mElements;

    public MdTag(final String tag)
    {
        mTag = tag;
        mElements = Lists.newArrayList();

        parseTag();
    }

    public List<MdTagElement> elements() { return mElements; }

    protected static final char DEL_TAG = '^';
    protected static final char MATCH_TAG = '0';

    public int baseLength()
    {
        // includes deleted bases
        return mElements.stream().mapToInt(x -> x.Length).sum();
    }

    public int sequenceLength()
    {
        return mElements.stream().filter(x -> x.Type != MdTagType.DEL).mapToInt(x -> x.Length).sum();
    }

    public boolean hasMismatches() { return mElements.stream().anyMatch(x -> x.Type != MdTagType.MATCH); }
    public int mismatchCount() { return (int)mElements.stream().filter(x -> x.Type != MdTagType.MATCH).count(); }

    public static final byte MATCH_BYTE = (byte)MdTagType.MATCH.ordinal();

    public byte[] extractSubSequence(final int seqIndexStart, final int seqIndexEnd, boolean reverse)
    {
        // build an array of values representing matches or mismatches for the sequence vs the ref genome from this alignment tag
        int seqIndex = 0;

        List<MdTagElement> elements;

        if(reverse)
        {
            elements = Lists.newArrayList(mElements);
            Collections.reverse(elements);
        }
        else
        {
            elements = mElements;
        }

        int length = seqIndexEnd - seqIndexStart + 1;
        byte[] matchArray = new byte[length];
        int index = 0;
        boolean followsDel = false;

        for(MdTagElement element : elements)
        {
            for(int i = 0; i < element.Length; ++i)
            {
                if(followsDel)
                {
                    followsDel = false;
                    ++seqIndex;
                    continue;
                }

                if(seqIndex > seqIndexEnd)
                    break;

                if(seqIndex >= seqIndexStart)
                {
                    matchArray[index++] = (byte)element.Type.ordinal();

                    if(element.Type == MdTagType.DEL)
                    {
                        followsDel = true;
                        break; // count this mismatch once only
                    }
                }

                ++seqIndex;
            }

            if(seqIndex > seqIndexEnd)
                break;
        }

        return matchArray;
    }

    private void parseTag()
    {
        int index = 0;
        int length = 0;
        MdTagType prevTag = null;

        while(index < mTag.length())
        {
            char c = mTag.charAt(index);

            MdTagType tag;

            if(c == DEL_TAG)
                tag = MdTagType.DEL;
            else if(isBase(c))
                tag = MdTagType.SNV;
            else
                tag = MdTagType.MATCH;

            if(prevTag == MdTagType.DEL && tag == MdTagType.SNV)
                tag = MdTagType.DEL; // since lists deleted bases

            if(prevTag != null && tag != prevTag)
            {
                if((prevTag == MdTagType.MATCH || prevTag == MdTagType.DEL) && length > 0)
                {
                    char prevChar = prevTag == MdTagType.MATCH ? MATCH_TAG : DEL_TAG;
                    mElements.add(new MdTagElement(prevTag, length, prevChar));
                }

                length = 0;

                if(tag == MdTagType.SNV)
                {
                    mElements.add(new MdTagElement(tag, 1, c));
                }
                else if(tag == MdTagType.MATCH)
                {
                    length = Integer.valueOf(String.valueOf(c));
                }
            }
            else if(prevTag == null || tag == prevTag)
            {
                if(tag == MdTagType.MATCH)
                    length = length * 10 + Integer.valueOf(String.valueOf(c));
                else if(tag == MdTagType.DEL)
                    ++length;
            }

            ++index;
            prevTag = tag;
        }

        if(prevTag == MdTagType.MATCH)
            mElements.add(new MdTagElement(prevTag, length, MATCH_TAG));
    }

    private static boolean isBase(char c)
    {
        return Nucleotides.baseIndex(c) >= 0;
    }

    public String toString()
    {
        return format("%s mismatches(snv=%d del=%d)",
                mTag, mElements.stream().filter(x -> x.isSnv()).count(), mElements.stream().filter(x -> x.isIndel()).count());
    }
}
