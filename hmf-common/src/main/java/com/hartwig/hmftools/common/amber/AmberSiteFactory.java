package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.VariantContext;

public final class AmberSiteFactory
{
    private static final Logger LOGGER = LogManager.getLogger(AmberSiteFactory.class);
    private static final String SNPCHECK = "SNPCHECK";

    public static AmberSite tumorSite(final TumorBAF baseDepth)
    {
        return ImmutableAmberSite.builder()
                .from(baseDepth)
                .snpCheck(false)
                .ref(baseDepth.ref())
                .alt(baseDepth.alt())
                .build();
    }

    public static AmberSite asSite(final BaseDepth baseDepth)
    {
        return ImmutableAmberSite.builder()
                .from(baseDepth)
                .snpCheck(false)
                .ref(baseDepth.ref().toString())
                .alt(baseDepth.alt().toString())
                .build();
    }

    public static ListMultimap<Chromosome, AmberSite> sites(final String filename) throws IOException
    {
        BufferedReader reader = createBufferedReader(filename);

        String header = reader.readLine();

        if(header.contains("fileformat=VCF"))
            return loadVcf(filename);
        else
            return loadFile(filename);
    }

    private static ListMultimap<Chromosome, AmberSite> loadVcf(final String vcfFile) throws IOException
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
            result.put(chromosome,
                    ImmutableAmberSite.builder()
                            .chromosome(variant.getContig())
                            .position(variant.getStart())
                            .ref(variant.getReference().getBaseString())
                            .alt(variant.getAlternateAllele(0).getBaseString())
                            .snpCheck(variant.hasAttribute(SNPCHECK))
                            .build());
        }

        LOGGER.info("loaded {} Amber germline sites", result.size());
        return result;
    }

    public static final String FLD_CHROMOSOME = "Chromosome";
    public static final String FLD_POSITION = "Position";
    public static final String FLD_REF = "Ref";
    public static final String FLD_ALT = "Alt";
    public static final String FLD_SNP_CHECK = "SnpCheck";

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

            result.put(chromosome,
                    ImmutableAmberSite.builder()
                            .chromosome(chrStr)
                            .position(Integer.parseInt(values[posIndex]))
                            .ref(values[refIndex])
                            .alt(values[altIndex])
                            .snpCheck(Boolean.parseBoolean(values[infoIndex]))
                            .build());
        }

        LOGGER.info("loaded {} Amber germline sites", result.size());
        return result;
    }
}
