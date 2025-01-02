package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Objects;

public class Chr7Genome extends SimpleTestGenome
{
    private static final int CHR7_EXTRACT_START = 140_000_000;
    private static final int CHR7_EXTRACT_LENGTH = 10_000_000;
    private final byte[] chr7Bytes;

    public Chr7Genome()
    {
        super();
        ClassLoader classLoader = getClass().getClassLoader();
        File file = new File(Objects.requireNonNull(classLoader.getResource("chr7_part.txt")).getFile());
        try
        {
            chr7Bytes = Files.readAllBytes(Path.of(file.getAbsolutePath()));
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    @Override
    public byte[] getBases(final String chromosome, final int posStart, final int posEnd)
    {
        if (posStart < CHR7_EXTRACT_START)
        {
            throw new IllegalArgumentException("posStart < " + CHR7_EXTRACT_START);
        }
        int length = posEnd - posStart + 1;
        if (length > CHR7_EXTRACT_LENGTH)
        {
            throw new IllegalArgumentException("length > " + CHR7_EXTRACT_LENGTH);
        }
        byte[] result = new byte[length];
        System.arraycopy(chr7Bytes, posStart - CHR7_EXTRACT_START, result, 0, length);
        return result;
    }

    @Override
    public String getBaseString(final String chromosome, final int posStart, final int posEnd)
    {
        return new String(getBases(chromosome, posStart, posEnd));
    }
}
