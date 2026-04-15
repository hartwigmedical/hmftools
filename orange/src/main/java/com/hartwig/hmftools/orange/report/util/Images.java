package com.hartwig.hmftools.orange.report.util;

import java.io.IOException;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.graphics.image.PDImageXObject;

public final class Images
{
    public static PDImageXObject load(final String path, final PDDocument document)
    {
        try
        {
            return PDImageXObject.createFromFile(path, document);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Could not read image from " + path, e);
        }
    }
}
