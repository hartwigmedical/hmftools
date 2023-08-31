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

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
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
            return loadVcf(filename);
        else
            return loadFile(filename);
    }

    public static ListMultimap<Chromosome, AmberSite> loadVcf(final String vcfFile) throws IOException
    {
        final ListMultimap<Chromosome, AmberSite> result = ArrayListMultimap.create();

        VcfFileReader reader = new VcfFileReader(vcfFile);

        if(!reader.fileValid())
            throw new IOException("invalid Amber sites file");

        for(VariantContext variant : reader.iterator())
        {
            if(variant.isFiltered())
                continue;

            if(!HumanChromosome.contains(variant.getContig()))
                continue;

            HumanChromosome chromosome = HumanChromosome.fromString(variant.getContig());

            result.put(chromosome, new AmberSite(
                    variant.getContig(), variant.getStart(), variant.getReference().getBaseString(),
                    variant.getAlternateAllele(0).getBaseString(), variant.hasAttribute(SNPCHECK)));
        }

        return result;
    }

    public static final String FLD_SNP_CHECK = "SnpCheck";

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(FLD_CHROMOSOME);
        sj.add(FLD_POSITION);
        sj.add(FLD_REF);
        sj.add(FLD_ALT);
        sj.add(FLD_SNP_CHECK);
        return sj.toString();
    }

    public static ListMultimap<Chromosome,AmberSite> loadFile(final String filename) throws IOException
    {
        final ListMultimap<Chromosome,AmberSite> result = ArrayListMultimap.create();

        BufferedReader reader = createBufferedReader(filename);

        String header = reader.readLine();

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
        int posIndex = fieldsIndexMap.get(FLD_POSITION);
        int refIndex = fieldsIndexMap.get(FLD_REF);
        int altIndex = fieldsIndexMap.get(FLD_ALT);
        int infoIndex = fieldsIndexMap.get(FLD_SNP_CHECK);

        String line = null;

        while((line = reader.readLine()) != null)
        {
            String[] values = line.split(TSV_DELIM, -1);

            String chrStr = values[chrIndex];

            if(!HumanChromosome.contains(chrStr))
                continue;

            HumanChromosome chromosome = HumanChromosome.fromString(chrStr);

            result.put(chromosome, new AmberSite(
                    chrStr, Integer.parseInt(values[posIndex]), values[refIndex], values[altIndex], Boolean.parseBoolean(values[infoIndex])));
        }

        LOGGER.info("loaded {} Amber germline sites from {}", result.size(), filename);
        return result;
    }
}
