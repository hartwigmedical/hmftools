package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

public class TinyGenome extends SimpleTestGenome
{
    private final Set<ChromosomeSnippet> data = new HashSet<>();

    public TinyGenome()
    {
        super();
        data.add( new ChromosomeSnippet("chr1", 10_000_000, 3_000_000, "chr1_part.txt"));
        data.add( new ChromosomeSnippet("chr3", 10_000_000, 3_000_000, "chr3_part.txt"));
        data.add( new ChromosomeSnippet("chr4", 105_000_000, 3_000_000, "chr4_part.txt"));
        data.add( new ChromosomeSnippet("chr5", 68_000_000, 1_000_000, "chr5_part_68.txt"));
        data.add( new ChromosomeSnippet("chr7", 140_000_000, 10_000_000, "chr7_part.txt"));
        data.add( new ChromosomeSnippet("chr7", 55_000_000, 3_000_000, "chr7_part_55.txt"));
    }

    @Override
    public byte[] getBases(final String chromosome, final int posStart, final int posEnd)
    {
        ChromosomeSnippet relevantSnippet = data.stream()
                .filter(snippet -> snippet.containsData(chromosome, posStart, posEnd))
                .findFirst().orElseThrow();
        return relevantSnippet.getBases(posStart, posEnd);
    }

    @Override
    public String getBaseString(final String chromosome, final int posStart, final int posEnd)
    {
        return new String(getBases(chromosome, posStart, posEnd));
    }
}

class ChromosomeSnippet
{
    private final int Start;
    private final int Length;
    private final int End;
    private final byte[] bytes;
    private final String Chromosome;

    ChromosomeSnippet(final String chromosome, final int start, final int length, final String dataFileName)
    {
        Start = start;
        Length = length;
        End = start + length;
        this.Chromosome = chromosome;
        ClassLoader classLoader = getClass().getClassLoader();
        File file = new File(Objects.requireNonNull(classLoader.getResource(dataFileName)).getFile());
        try
        {
            bytes = Files.readAllBytes(Path.of(file.getAbsolutePath()));
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    boolean containsData(String chromosome, int start, int end)
    {
        if (!chromosome.equals(Chromosome))
        {
            return false;
        }
        return start >= Start && end <= End;
    }

    public byte[] getBases(final int posStart, final int posEnd)
    {
        if (posStart < Start)
        {
            throw new IllegalArgumentException("posStart < " + Start);
        }
        int length = posEnd - posStart + 1;
        if (length > Length)
        {
            throw new IllegalArgumentException("length > " + Length);
        }
        byte[] result = new byte[length];
        System.arraycopy(bytes, posStart - Start, result, 0, length);
        return result;
    }
}
