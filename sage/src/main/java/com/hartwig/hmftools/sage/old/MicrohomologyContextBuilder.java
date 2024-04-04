package com.hartwig.hmftools.sage.old;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;

public final class MicrohomologyContextBuilder
{
    public static String microhomologyAtInsert(int position, final String refSequence, final String alt)
    {
        if(refSequence.length() < position)
        {
            SG_LOGGER.warn("Attempt to determine microhomology outside of sequence length");
            return Strings.EMPTY;
        }

        if(alt.contains(","))
        {
            return Strings.EMPTY;
        }

        final String readSequence = refSequence.substring(0, position) + alt + refSequence.substring(position + 1);
        return microhomologyAtInsert(position, alt.length(), readSequence.getBytes()).toString();
    }

    public static String microhomologyAtDelete(int position, final String refSequence, final String ref)
    {
        if(refSequence.length() < position + ref.length())
        {
            SG_LOGGER.warn("Attempt to determine microhomology outside of sequence length");
            return Strings.EMPTY;
        }

        return microhomologyAtDelete(position, ref.length(), refSequence.getBytes()).toString();
    }

    public static MicrohomologyContext microhomologyAtDeleteFromReadSequence(int position, final String ref,
            final byte[] readSequence)
    {
        return microhomologyAtDelete(position, ref.length(), reconstructDeletedSequence(position, readSequence, ref));
    }

    public static MicrohomologyContext microhomologyAtInsert(int position, int altLength, final byte[] readSequence)
    {
        return leftAlignedMicrohomology(position, altLength, readSequence);
    }

    public static MicrohomologyContext microhomologyAtDelete(int position, int refLength, final byte[] sequence)
    {
        return leftAlignedMicrohomology(position, refLength, sequence);
    }

    private static MicrohomologyContext leftAlignedMicrohomology(int position, int altLength, final byte[] sequence)
    {
        // Left align
        if(position > 0 && position + altLength - 1 < sequence.length)
        {
            byte ref = sequence[position];
            byte alt = sequence[position + altLength - 1];
            if(ref == alt)
            {
                return microhomologyAtInsert(position - 1, altLength, sequence);
            }
        }

        int length = 0;
        for(int i = 0; i < altLength - 1; i++)
        {
            if(position + i + altLength < sequence.length && sequence[position + i + 1] == sequence[position + i + altLength])
            {
                length++;
            }
            else
            {
                break;
            }
        }

        return new MicrohomologyContext(position, sequence, length);
    }

    public static MicrohomologyContext expandMicrohomologyRepeats(final MicrohomologyContext result)
    {
        if(result.length() == 0)
        {
            return result;
        }

        final byte[] readSequence = result.readSequence();
        int length = result.length();
        int mhIndex = 0;
        for(int i = result.homologyIndex() + length; i < readSequence.length; i++)
        {
            if(readSequence[i] == readSequence[result.homologyIndex() + mhIndex])
            {
                length++;
                mhIndex++;
                if(mhIndex >= result.length())
                {
                    mhIndex = 0;
                }
            }
            else
            {
                break;
            }
        }

        if(result.length() != length)
        {
            return new MicrohomologyContext(result.position(), readSequence, length);
        }

        return result;
    }

    @VisibleForTesting
    public static byte[] reconstructDeletedSequence(int position, final byte[] readSequence, final String ref)
    {
        final byte[] refBytes = ref.getBytes();
        final byte[] completeSequence = new byte[readSequence.length + ref.length() - 1];

        int tailLength = readSequence.length - position - 1;

        System.arraycopy(readSequence, 0, completeSequence, 0, position + 1);
        System.arraycopy(refBytes, 1, completeSequence, position + 1, ref.length() - 1);
        System.arraycopy(readSequence, position + 1, completeSequence, completeSequence.length - tailLength, tailLength);

        return completeSequence;
    }
}
