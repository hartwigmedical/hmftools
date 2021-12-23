package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.common.utils.FileWriterUtils;

public class FileSources
{
    public final String Source;
    public final String Linx;
    public final String Somatic;
    public final String Purple;

    public static final String SAMPLE_DIR = "sample_dir";
    public static final String LINX_DIR = "linx_dir";
    public static final String PURPLE_DIR = "purple_dir";
    public static final String SOMATIC_DIR = "somatic_dir";

    public FileSources(final String source, final String linx, final String somatic, final String purple)
    {
        Source = source;
        Linx = linx;
        Somatic = somatic;
        Purple = purple;
    }

    public static FileSources fromConfig(final String fileSourceStr)
    {
        String[] values = fileSourceStr.split(ITEM_DELIM, -1);

        if(values.length < 2)
        {
            CMP_LOGGER.error("invalid file source config entry", fileSourceStr);
            return null;
        }

        String source = values[0];
        String sampleDir = "";
        String linxDir = "";
        String purpleDir = "";
        String somaticDir = "";

        int itemIndex = 1;

        if(values[itemIndex].startsWith(SAMPLE_DIR))
        {
            sampleDir = checkAddDirSeparator(values[1].split(SUB_ITEM_DELIM)[1]);
            ++itemIndex;
        }

        for(int i = itemIndex; i < values.length; ++i)
        {
            String[] itemStr = values[i].split(SUB_ITEM_DELIM);

            if(itemStr.length != 2)
                return null;

            if(itemStr[0].equals(LINX_DIR))
            {
                linxDir = getDirectory(sampleDir, itemStr[1]);
            }
            else if(itemStr[0].equals(PURPLE_DIR))
            {
                purpleDir = getDirectory(sampleDir, itemStr[1]);
            }
            else if(itemStr[0].equals(SOMATIC_DIR))
            {
                somaticDir = getDirectory(sampleDir, itemStr[1]);
            }
        }

        if(linxDir.isEmpty())
            linxDir = sampleDir;

        if(purpleDir.isEmpty())
            purpleDir = sampleDir;

        if(somaticDir.isEmpty())
            somaticDir = sampleDir;

        return new FileSources(source, linxDir, purpleDir, somaticDir);
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
