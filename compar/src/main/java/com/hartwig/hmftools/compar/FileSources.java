package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.SUB_ITEM_DELIM;

import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;

public class FileSources
{
    public final String Source;
    public final String Linx;
    public final String LinxGermline;
    public final String Purple;
    public final String Cuppa;
    public final String Lilac;
    public final String Chord;

    public static final String SAMPLE_DIR = "sample_dir";
    public static final String LINX_DIR = "linx_dir";
    public static final String LINX_GERMLINE_DIR = "linx_germline_dir";
    public static final String PURPLE_DIR = "purple_dir";
    public static final String CUPPA_DIR = "cuppa_dir";
    public static final String LILAC_DIR = "lilac_dir";
    public static final String CHORD_DIR = "chord_dir";

    public FileSources(final String source, final String linx, final String purple, final String linxGermline, final String cuppa,
            final String lilac, final String chord)
    {
        Source = source;
        Linx = linx;
        LinxGermline = linxGermline;
        Purple = purple;
        Cuppa = cuppa;
        Lilac = lilac;
        Chord = chord;
    }

    public static FileSources sampleInstance(final FileSources fileSources, final String sampleId)
    {
        return new FileSources(
                fileSources.Source,
                fileSources.Linx.replaceAll("\\*", sampleId),
                fileSources.Purple.replaceAll("\\*", sampleId),
                fileSources.LinxGermline.replaceAll("\\*", sampleId),
                fileSources.Cuppa.replaceAll("\\*", sampleId),
                fileSources.Lilac.replaceAll("\\*", sampleId),
                fileSources.Chord.replaceAll("\\*", sampleId));
    }

    public static FileSources fromConfig(final String sourceName, final String fileSourceStr)
    {
        String[] values = fileSourceStr.split(ITEM_DELIM, -1);

        String sampleDir = "";

        int itemIndex = 0;

        if(values[itemIndex].startsWith(SAMPLE_DIR))
        {
            sampleDir = checkAddDirSeparator(values[itemIndex].split(SUB_ITEM_DELIM)[1]);
            ++itemIndex;
        }

        // by default use the sample root directory and then default pipeline directory names per tool and mode
        String linxDir = getDirectory(sampleDir, PipelineToolDirectories.LINX_SOMATIC_DIR);
        String linxGermlineDir = getDirectory(sampleDir, PipelineToolDirectories.LINX_GERMLINE_DIR);
        String purpleDir = getDirectory(sampleDir, PipelineToolDirectories.PURPLE_DIR);
        String cuppaDir = getDirectory(sampleDir, PipelineToolDirectories.CUPPA_DIR);
        String lilacDir = getDirectory(sampleDir, PipelineToolDirectories.LILAC_DIR);
        String chordDir = getDirectory(sampleDir, PipelineToolDirectories.CHORD_DIR);

        for(int i = itemIndex; i < values.length; ++i)
        {
            String[] itemStr = values[i].split(SUB_ITEM_DELIM);

            if(itemStr.length != 2)
                return null;

            String type = itemStr[0];
            String value = itemStr[1];

            if(type.equals(LINX_DIR))
            {
                linxDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(LINX_GERMLINE_DIR))
            {
                linxGermlineDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(PURPLE_DIR))
            {
                purpleDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(LILAC_DIR))
            {
                lilacDir = value;
            }
            else if(type.equals(CHORD_DIR))
            {
                chordDir = value;
            }
            else if(type.equals(CUPPA_DIR))
            {
                cuppaDir = value;
            }
        }

        return new FileSources(sourceName, linxDir, purpleDir, linxGermlineDir, cuppaDir, lilacDir, chordDir);
    }

    private static String getDirectory(final String sampleDir, final String typeDir)
    {
        String directory = "";

        if(sampleDir.isEmpty())
            directory = typeDir;
        else if(typeDir.isEmpty())
            directory = sampleDir;
        else
            directory = String.format("%s%s", sampleDir, typeDir);

        return checkAddDirSeparator(directory);
    }
}
