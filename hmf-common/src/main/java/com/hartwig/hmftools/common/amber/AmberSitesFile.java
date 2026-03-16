package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.VariantContext;

public final class AmberSitesFile
{
    private static final Logger LOGGER = LogManager.getLogger(AmberSitesFile.class);
    private static final String SNPCHECK = "SNPCHECK";

    public static ListMultimap<Chromosome,AmberSite> sites(final String filename) throws IOException
    {
        BufferedReader reader = createBufferedReader(filename);

        String header = reader.readLine();

        if(header.contains("fileformat=VCF"))
        {
            return loadVcf(filename);
        }
        else
        {
            return loadFile(filename);
        }
    }

    public static ListMultimap<Chromosome, AmberSite> loadVcf(final String vcfFile) throws IOException
    {
        final ListMultimap<Chromosome, AmberSite> result = ArrayListMultimap.create();

        VcfFileReader reader = new VcfFileReader(vcfFile);

        if(!reader.fileValid())
        {
            throw new IOException("invalid Amber sites file");
        }

        for(VariantContext variant : reader.iterator())
        {
            if(variant.isFiltered())
            {
                continue;
            }

            if(!HumanChromosome.contains(variant.getContig()))
            {
                continue;
            }

            HumanChromosome chromosome = HumanChromosome.fromString(variant.getContig());

            result.put(chromosome, new AmberSite(
                    variant.getContig(), variant.getStart(), variant.getReference().getBaseString(),
                    variant.getAlternateAllele(0).getBaseString(), variant.hasAttribute(SNPCHECK)));
        }
        return result;
    }

    public static final String FLD_SNP_CHECK = "SnpCheck";
    public static final String FLD_FREQUENCY = "Frequency";

    public static List<String> columns()
    {
        return List.of(FLD_CHROMOSOME, FLD_POSITION, FLD_REF, FLD_ALT, FLD_SNP_CHECK, FLD_FREQUENCY);
    }

    public static String header()
    {
        return String.join(TSV_DELIM, columns());
    }

    public static ListMultimap<Chromosome, AmberSite> loadFile(final String filename) throws IOException
    {
        final ListMultimap<Chromosome, AmberSite> result = ArrayListMultimap.create();

        BufferedReader reader = createBufferedReader(filename);

        String header = reader.readLine();

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
        int posIndex = fieldsIndexMap.get(FLD_POSITION);
        int refIndex = fieldsIndexMap.get(FLD_REF);
        int altIndex = fieldsIndexMap.get(FLD_ALT);
        int infoIndex = fieldsIndexMap.get(FLD_SNP_CHECK);

        String line;

        while((line = reader.readLine()) != null)
        {
            String[] values = line.split(TSV_DELIM, -1);

            String chrStr = values[chrIndex];

            if(!HumanChromosome.contains(chrStr))
            {
                continue;
            }

            HumanChromosome chromosome = HumanChromosome.fromString(chrStr);

            final int position = Integer.parseInt(values[posIndex]);
            final boolean snpCheck = Boolean.parseBoolean(values[infoIndex]);
            final String ref = values[refIndex];
            final String alt = values[altIndex];
            if(values.length == 5)
            {
                result.put(chromosome, new AmberSite(chrStr, position, ref, alt, snpCheck));
            }
            else
            {
                double gnomadFrequency = Double.parseDouble(values[5]);
                result.put(chromosome, new AmberSite(chrStr, position, ref, alt, snpCheck, gnomadFrequency));
            }
        }
        LOGGER.info("loaded {} Amber germline sites from {}", result.size(), filename);
        return result;
    }

    public static void writeData(final ListMultimap<Chromosome, AmberSite> data, final String absolutePath)
    {
        SortedSet<HumanChromosome> chromosomesInOrder = new TreeSet<>();
        for(Chromosome chromosome : data.asMap().keySet())
        {
            chromosomesInOrder.add((HumanChromosome) chromosome);
        }
        DelimFileWriter<AmberSite> writer = new DelimFileWriter<>(absolutePath, columns(), new AmberSiteEncoder());
        for(Chromosome chromosome : chromosomesInOrder)
        {
            List<AmberSite> sites = data.get(chromosome);
            for(AmberSite site : sites)
            {
                writer.writeRow(site);
            }
        }
        writer.close();
    }

    static class AmberSiteEncoder implements BiConsumer<AmberSite, DelimFileWriter.Row>
    {
        @Override
        public void accept(final AmberSite amberSite, final DelimFileWriter.Row row)
        {
            row.set(FLD_CHROMOSOME, amberSite.Chromosome);
            row.set(FLD_POSITION, amberSite.Position);
            row.set(FLD_REF, amberSite.Ref);
            row.set(FLD_ALT, amberSite.Alt);
            row.set(FLD_SNP_CHECK, amberSite.snpCheck());
            row.set(FLD_FREQUENCY, amberSite.GnomadFrequency);
        }
    }
}
