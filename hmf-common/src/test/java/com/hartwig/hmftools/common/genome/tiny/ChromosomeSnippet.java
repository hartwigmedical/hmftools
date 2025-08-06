package com.hartwig.hmftools.common.genome.tiny;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Objects;

public class ChromosomeSnippet
{
    private final int Start;
    private final int Length;
    private final int End;
    private final byte[] bytes;
    private final String Chromosome;

    public ChromosomeSnippet(final String chromosome, final int start, final int length, final Path dataFileName)
    {
        Start = start;
        Length = length;
        End = start + length;
        this.Chromosome = chromosome;
        try
        {
            bytes = Files.readAllBytes(dataFileName);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    boolean containsData(String chromosome, int start, int end)
    {
        if(!chromosome.equals(Chromosome))
        {
            return false;
        }
        return start >= Start && end <= End;
    }

    public byte[] getBases(final int posStart, final int posEnd)
    {
        if(posStart < Start)
        {
            throw new IllegalArgumentException("posStart < " + Start);
        }
        int length = posEnd - posStart + 1;
        if(length > Length)
        {
            throw new IllegalArgumentException("length > " + Length);
        }
        byte[] result = new byte[length];
        System.arraycopy(bytes, posStart - Start, result, 0, length);
        return result;
    }
}
