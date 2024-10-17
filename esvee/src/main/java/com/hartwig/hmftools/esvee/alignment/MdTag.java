package com.hartwig.hmftools.esvee.alignment;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;

public class MdTag
{
    private final String mTag;
    private final List<MdTagElement> mElements;
    private final int mBaseLength;

    public MdTag(final String tag)
    {
        mTag = tag;
        mElements = Lists.newArrayList();

        parseTag();

        mBaseLength = mElements.stream().mapToInt(x -> x.Length).sum();
    }

    public List<MdTagElement> elements() { return mElements; }

    protected static final char DEL_TAG = '^';
    protected static final char MATCH_TAG = '0';

    public int baseLength() { return mBaseLength; }

    public boolean hasMismatches() { return mElements.stream().anyMatch(x -> x.Type != MdTagType.MATCH); }

    public byte[] extractSubSequence(final int seqIndexStart, final int seqIndexEnd, boolean reverse)
    {
        // establish the sequence start of the MD tag
        StringBuilder sb = new StringBuilder();

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

        for(MdTagElement element : elements)
        {
            for(int i = 0; i < element.Length; ++i)
            {
                if(element.Type != MdTagType.DEL)
                {
                    if(seqIndex >= seqIndexStart)
                        sb.append(element.Base);

                    ++seqIndex;
                }

                if(seqIndex > seqIndexEnd)
                    break;
            }

            if(seqIndex > seqIndexEnd)
                break;
        }

        return sb.toString().getBytes();
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
}
