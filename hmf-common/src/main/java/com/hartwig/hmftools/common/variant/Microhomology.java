package com.hartwig.hmftools.common.variant;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class Microhomology
{
    private static final Logger LOGGER = LogManager.getLogger(Microhomology.class);

    @NotNull
    public static String microhomologyAtInsert(int position, @NotNull final String refSequence, @NotNull final String alt)
    {
        if(refSequence.length() < position)
        {
            LOGGER.warn("Attempt to determine microhomology outside of sequence length");
            return Strings.EMPTY;
        }

        if(alt.contains(","))
        {
            return Strings.EMPTY;
        }

        final String readSequence = refSequence.substring(0, position) + alt + refSequence.substring(position + 1);
        return microhomologyAtInsert(position, alt.length(), readSequence.getBytes()).toString();
    }

    @NotNull
    public static String microhomologyAtDelete(int position, @NotNull final String refSequence, @NotNull final String ref)
    {
        if(refSequence.length() < position + ref.length())
        {
            LOGGER.warn("Attempt to determine microhomology outside of sequence length");
            return Strings.EMPTY;
        }

        return microhomologyAtDelete(position, ref.length(), refSequence.getBytes()).toString();
    }

    @NotNull
    public static MicrohomologyContext microhomologyAtDeleteFromReadSequence(int position, @NotNull final String ref,
            @NotNull final byte[] readSequence)
    {
        return microhomologyAtDelete(position, ref.length(), reconstructDeletedSequence(position, readSequence, ref));
    }

    @NotNull
    public static MicrohomologyContext microhomologyAtInsert(int position, int altLength, @NotNull final byte[] readSequence)
    {
        return leftAlignedMicrohomology(position, altLength, readSequence);
    }

    @NotNull
    public static MicrohomologyContext microhomologyAtDelete(int position, int refLength, @NotNull final byte[] sequence)
    {
        return leftAlignedMicrohomology(position, refLength, sequence);
    }

    private static MicrohomologyContext leftAlignedMicrohomology(int position, int altLength, @NotNull final byte[] sequence)
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

    @NotNull
    public static MicrohomologyContext expandMicrohomologyRepeats(@NotNull final MicrohomologyContext result)
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
    static byte[] reconstructDeletedSequence(int position, @NotNull final byte[] readSequence, @NotNull final String ref)
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
