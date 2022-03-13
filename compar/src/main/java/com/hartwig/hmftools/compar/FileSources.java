package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

public class FileSources
{
    public final String Source;
    public final String Linx;
    public final String Purple;
    public final String SomaticVcf;
    public final String GermlineVcf;

    public static final String SAMPLE_DIR = "sample_dir";
    public static final String LINX_DIR = "linx_dir";
    public static final String PURPLE_DIR = "purple_dir";

    // specify Sage VCFs instead of the default Purple VCFs
    public static final String SOMATIC_VCF = "somatic_vcf";
    public static final String GERMLINE_VCF = "germline_vcf";

    public FileSources(final String source, final String linx, final String purple, final String somaticVcf, final String germlineVcf)
    {
        Source = source;
        Linx = linx;
        Purple = purple;
        SomaticVcf = somaticVcf;
        GermlineVcf = germlineVcf;
    }

    public static FileSources sampleInstance(final FileSources fileSources, final String sampleId)
    {
        if(!fileSources.Linx.contains("*") && !fileSources.Purple.contains("*") && !fileSources.SomaticVcf.contains("*")
        && !fileSources.GermlineVcf.contains("*"))
        {
            return fileSources;
        }

        return new FileSources(
                fileSources.Source,
                fileSources.Linx.replaceAll("\\*", sampleId),
                fileSources.Purple.replaceAll("\\*", sampleId),
                fileSources.SomaticVcf.replaceAll("\\*", sampleId),
                fileSources.GermlineVcf.replaceAll("\\*", sampleId));
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
        String somaticVcf = "";
        String germlineVcf = "";

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

            String type = itemStr[0];
            String value = itemStr[1];

            if(type.equals(LINX_DIR))
            {
                linxDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(PURPLE_DIR))
            {
                purpleDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(SOMATIC_VCF))
            {
                somaticVcf = value;
            }
            else if(type.equals(GERMLINE_VCF))
            {
                germlineVcf = value;
            }
        }

        if(linxDir.isEmpty())
            linxDir = sampleDir;

        if(purpleDir.isEmpty())
            purpleDir = sampleDir;

        return new FileSources(source, linxDir, purpleDir, somaticVcf, germlineVcf);
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
