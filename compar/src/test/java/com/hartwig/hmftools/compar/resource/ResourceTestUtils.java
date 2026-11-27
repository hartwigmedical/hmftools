package com.hartwig.hmftools.compar.resource;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public final class ResourceTestUtils
{
    public static File writeFile(final File dir, final String relativePath, final String... lines) throws IOException
    {
        File file = new File(dir, relativePath);
        file.getParentFile().mkdirs();

        try(BufferedWriter writer = createBufferedWriter(file.getAbsolutePath()))
        {
            for(String line : lines)
            {
                writer.write(line);
                writer.newLine();
            }
        }

        return file;
    }

    public static String tsv(final String... values)
    {
        return String.join(TSV_DELIM, values);
    }
}
