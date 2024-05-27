package com.hartwig.hmftools.bamtools.compare;

import java.io.File;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CompareUtils
{
    public static SamReaderFactory makeSamReaderFactory(CompareConfig config)
    {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT);
        if(config.RefGenomeFile != null && !config.RefGenomeFile.isEmpty())
        {
            readerFactory.referenceSequence(new File(config.RefGenomeFile));
        }
        return readerFactory;
    }
}
