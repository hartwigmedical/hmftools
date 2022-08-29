package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

public class FileSources
{
    public final String Source;
    public final String Linx;
    public final String LinxGermline;
    public final String Purple;
    public final String SomaticVcf;
    public final String GermlineVcf;
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

    // specify Sage VCFs instead of the default Purple VCFs
    public static final String SOMATIC_VCF = "somatic_vcf";
    public static final String GERMLINE_VCF = "germline_vcf";

    public FileSources(final String source, final String linx, final String purple, final String somaticVcf, final String germlineVcf,
            final String linxGermline, final String cuppa, final String lilac, final String chord)
    {
        Source = source;
        Linx = linx;
        LinxGermline = linxGermline;
        Purple = purple;
        SomaticVcf = somaticVcf;
        GermlineVcf = germlineVcf;
        Cuppa = cuppa;
        Lilac = lilac;
        Chord = chord;
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
                fileSources.GermlineVcf.replaceAll("\\*", sampleId),
                fileSources.LinxGermline.replaceAll("\\*", sampleId),
                fileSources.Cuppa.replaceAll("\\*", sampleId),
                fileSources.Lilac.replaceAll("\\*", sampleId),
                fileSources.Chord.replaceAll("\\*", sampleId));
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
        String linxGermlineDir = "";
        String purpleDir = "";
        String somaticVcf = "";
        String germlineVcf = "";
        String cuppaDir = "";
        String lilacDir = "";
        String chordDir = "";

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
            else if(type.equals(LINX_GERMLINE_DIR))
            {
                linxGermlineDir = getDirectory(sampleDir, value);
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

        if(linxDir.isEmpty())
            linxDir = sampleDir;

        if(linxGermlineDir.isEmpty())
            linxGermlineDir = sampleDir;

        if(purpleDir.isEmpty())
            purpleDir = sampleDir;

        if(chordDir.isEmpty())
            chordDir = sampleDir;

        if(cuppaDir.isEmpty())
            cuppaDir = sampleDir;

        if(lilacDir.isEmpty())
            lilacDir = sampleDir;

        return new FileSources(source, linxDir, purpleDir, somaticVcf, germlineVcf, linxGermlineDir, cuppaDir, lilacDir, chordDir);
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
