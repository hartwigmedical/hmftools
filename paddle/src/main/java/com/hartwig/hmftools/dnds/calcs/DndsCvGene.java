package com.hartwig.hmftools.dnds.calcs;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class DndsCvGene
{
    public final String Gene;
    public final int Synonymous;
    public final DndsCv Missense;
    public final DndsCv Nonsense;
    public final DndsCv Splice;
    public final DndsCv Indels;

    public DndsCvGene(
            final String gene, final int synonymous, final DndsCv missense, final DndsCv nonsense, final DndsCv splice, final DndsCv indels)
    {
        Gene = gene;
        Synonymous = synonymous;
        Missense = missense;
        Nonsense = nonsense;
        Splice = splice;
        Indels = indels;
    }

    public static final String REF_CDSC_FILE = "HmfRefCDSCv.tsv";

    public static String geneDndsCvFilename(final String sourceDir)
    {
        return sourceDir + REF_CDSC_FILE;
    }

    public static Map<String,DndsCvGene> load(final String filename)
    {
        Map<String,DndsCvGene> geneData = Maps.newHashMap();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String gene = values[0];

                geneData.put(gene, new DndsCvGene(
                        gene, Integer.parseInt(values[1]),
                        new DndsCv(Integer.parseInt(values[2]), Double.parseDouble(values[6])),
                        new DndsCv(Integer.parseInt(values[3]), Double.parseDouble(values[7])),
                        new DndsCv(Integer.parseInt(values[4]), Double.parseDouble(values[8])),
                        new DndsCv(Integer.parseInt(values[5]), Double.parseDouble(values[9]))));
            }
        }
        catch(IOException e)
        {
            DN_LOGGER.error("failed to read DNDS CvGene file: {}", e.toString());
            return null;
        }

        return geneData;
    }
}
