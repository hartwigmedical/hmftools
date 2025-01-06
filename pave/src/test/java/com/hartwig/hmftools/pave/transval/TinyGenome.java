package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public class TinyGenome extends SimpleTestGenome
{
    private final Map<String, ChromosomeSnippet> data = new HashMap<>();

    public TinyGenome()
    {
        super();
        data.put("chr1", new ChromosomeSnippet(10_000_000, 3_000_000, "chr1_part.txt"));
        data.put("chr7", new ChromosomeSnippet(140_000_000, 10_000_000, "chr7_part.txt"));
    }

    @Override
    public byte[] getBases(final String chromosome, final int posStart, final int posEnd)
    {
        if (!data.containsKey(chromosome)) {
            throw new IllegalArgumentException("Unknown chromosome: " + chromosome);
        }
        return data.get(chromosome).getBases(posStart, posEnd);
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
    private final byte[] bytes;

    ChromosomeSnippet(final int start, final int length, final String dataFileName)
    {
        Start = start;
        Length = length;
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
