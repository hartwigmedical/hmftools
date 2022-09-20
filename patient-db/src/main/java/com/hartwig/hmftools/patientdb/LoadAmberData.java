package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberMapping;
import com.hartwig.hmftools.common.amber.AmberMappingFactory;
import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.common.amber.AmberPatientFactory;
import com.hartwig.hmftools.common.amber.AmberSample;
import com.hartwig.hmftools.common.amber.AmberSampleFactory;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSiteFactory;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFFileReader;

public class LoadAmberData
{
    private static final Logger LOGGER = LogManager.getLogger(LoadAmberData.class);

    private static final int DEFAULT_MIN_DEPTH = 10;
    private static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    private static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    private static final String SAMPLE = "sample";
    private static final String AMBER_SNP_VCF = "amber_snp_vcf";
    private static final String SNPCHECK_VCF = "snpcheck_vcf";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String amberSnpPath = cmd.getOptionValue(AMBER_SNP_VCF);
        String mappingLoci = cmd.getOptionValue(SNPCHECK_VCF);

        LOGGER.info("Loading mapping loci from {}", mappingLoci);
        final ListMultimap<Chromosome, AmberSite> mappingSites = AmberSiteFactory.sites(mappingLoci);

        final GenomePositionSelector<AmberSite> selector = GenomePositionSelectorFactory.create(mappingSites);

        try(final DatabaseAccess dbAccess = databaseAccess(cmd);
                final VCFFileReader fileReader = new VCFFileReader(new File(amberSnpPath), false))
        {
            LOGGER.info("Loading vcf snp data from {}", amberSnpPath);
            final List<BaseDepth> baseDepths = fileReader.iterator()
                    .stream()
                    .map(BaseDepthFactory::fromVariantContext)
                    .filter(x -> selector.select(x).isPresent())
                    .collect(Collectors.toList());

            final AmberSampleFactory amberSampleFactory =
                    new AmberSampleFactory(DEFAULT_MIN_DEPTH, DEFAULT_MIN_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);
            final AmberSample sample = amberSampleFactory.fromBaseDepth(tumorSample, baseDepths);

            processSample(sample, dbAccess);
        }

        LOGGER.info("Complete");
    }

    public static void processSample(final AmberSample sample, final DatabaseAccess dbAccess)
    {
        LOGGER.info("Comparing with existing samples");
        final List<AmberSample> allSamples = dbAccess.readAmberSamples();
        final List<AmberMapping> sampleMappings = Lists.newArrayList();
        for(AmberSample other : allSamples)
        {
            if(!other.sampleId().equals(sample.sampleId()))
            {
                final AmberMapping mapping = AmberMappingFactory.create(sample, other);
                if(mapping.likelihood() > 0.8)
                {
                    sampleMappings.add(mapping);
                }
            }
        }

        LOGGER.info("Sample {} matched with {} other samples", sample.sampleId(), sampleMappings.size());
        final List<AmberPatient> existingPatients = dbAccess.readAmberPatients();
        final AmberPatientFactory amberPatientFactory = new AmberPatientFactory(existingPatients, sampleMappings);
        final AmberPatient patient = amberPatientFactory.createPatient(sample);

        LOGGER.info("Writing data");
        dbAccess.writeAmberSample(sample);
        dbAccess.writeAmberMapping(sample.sampleId(), sampleMappings);
        dbAccess.writeAmberPatients(Collections.singletonList(patient));
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample");
        options.addOption(AMBER_SNP_VCF, true, "Path to the amber snp vcf");
        options.addOption(SNPCHECK_VCF, true, "Path to the downsampled snp check vcf");
        addDatabaseCmdLineArgs(options);
        return options;
    }
}
